"""One local, reproducible entry point for genotype evidence annotation."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import UTC, datetime
import hashlib
import json
from pathlib import Path
import platform
import sys
from typing import Any

from .analysis import build_analysis_report, print_console_summary, write_json_report, write_markdown_report
from .evidence import load_snapshot
from .utils import load_genome_data

VERSION = "0.2.0"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _code_hash() -> str:
    digest = hashlib.sha256()
    package = Path(__file__).parent
    for path in sorted(package.rglob('*.py')):
        if any(part in {'tests', '__pycache__', 'evidence_cache'} for part in path.parts):
            continue
        digest.update(path.relative_to(package).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _paths(args) -> tuple[Path, Path, Path]:
    output_dir = Path(args.output_dir).expanduser()
    prefix = args.output_prefix or f"gene_report_{datetime.now(UTC):%Y%m%d_%H%M%S_%f}"
    if Path(prefix).name != prefix or prefix in {'.', '..'}:
        raise ValueError('--output-prefix must be a simple filename, without directories.')
    return (Path(args.json_output).expanduser() if args.json_output else output_dir / f'{prefix}.json',
            Path(args.markdown_output).expanduser() if args.markdown_output else output_dir / f'{prefix}.md',
            output_dir / f'{prefix}.manifest.json')


def _get_clinical(args, genome):
    from .evidence.clinical import load_clinvar_manifest, match_clinvar, refresh_clinvar_database, is_clinvar_stale
    assembly = genome.assembly
    # Validate a pinned artifact even when input provenance prevents matching.
    manifest = None
    if args.clinvar_manifest:
        manifest = load_clinvar_manifest(args.clinvar_manifest,
            assembly=assembly if assembly in {'GRCh37', 'GRCh38'} else None)
    if assembly not in {'GRCh37', 'GRCh38'}:
        return {'status': 'not_analyzed', 'reason': 'Genome build is unknown. Confirm the build from your provider and use --assembly.', 'findings': []}
    if genome.strand != 'plus':
        return {'status': 'not_analyzed', 'reason': 'Plus-strand orientation is required. Confirm it with your provider and use --strand plus.', 'findings': []}
    if manifest is None:
        manifest = load_clinvar_manifest(assembly=assembly, cache_dir=args.cache_dir)
    if args.clinvar_manifest and manifest is None:
        raise FileNotFoundError(f'ClinVar manifest not found: {args.clinvar_manifest}')
    if not args.clinvar_manifest and (manifest is None or is_clinvar_stale(manifest, args.stale_after_days)):
        if args.auto_update and not args.offline:
            print(f'Updating the public ClinVar database ({assembly}); your genotype file stays local...')
            manifest = refresh_clinvar_database(assembly=assembly, cache_dir=args.cache_dir, show_progress=True)
        elif manifest is None:
            if args.snapshot_path:
                return {'status': 'not_analyzed', 'reason': 'Only historical marker evidence supplied; no ClinVar database available.', 'findings': []}
            raise FileNotFoundError('No ClinVar database is cached. Run update-evidence first, or pass --clinvar-manifest.')
    clinical = match_clinvar(genome, manifest_path=manifest['manifest_path'], assembly=assembly)
    clinical['cache_stale'] = is_clinvar_stale(manifest, args.stale_after_days)
    clinical['pinned'] = bool(args.clinvar_manifest)
    return clinical


def _get_gene_validity(args, clinical):
    from .evidence.gene_validity import load_gene_validity, refresh_gene_validity, attach_gene_validity, is_gene_validity_stale
    snapshot = load_gene_validity(path=args.clingen_snapshot, cache_dir=args.cache_dir)
    if args.clingen_snapshot and snapshot is None:
        raise FileNotFoundError(f'ClinGen snapshot not found: {args.clingen_snapshot}')
    if not args.clingen_snapshot and (snapshot is None or is_gene_validity_stale(snapshot, args.stale_after_days)):
        if args.auto_update and not args.offline:
            print('Updating the public ClinGen gene-disease validity table...')
            snapshot = refresh_gene_validity(cache_dir=args.cache_dir)
    if snapshot is None:
        return {'status': 'not_analyzed', 'reason': 'No ClinGen snapshot available in offline mode.'}
    attach_gene_validity(clinical, snapshot)
    return {'status': 'annotated', 'metadata': {**snapshot.get('metadata', {}), 'snapshot_sha256': snapshot.get('snapshot_sha256')},
            'cache_stale': is_gene_validity_stale(snapshot, args.stale_after_days),
            'scope': 'Gene-disease validity is contextual evidence, not proof that a matched variant causes disease.'}


def cmd_analyze(args) -> int:
    input_path = Path(args.genome_file).expanduser().resolve()
    outputs = _paths(args)
    if len({p.resolve() for p in outputs}) != 3 or input_path in {p.resolve() for p in outputs}:
        raise ValueError('Input and the three report outputs must have different paths.')
    if any(p.exists() for p in outputs):
        raise FileExistsError('A report output already exists. Choose a different --output-prefix or output path.')
    print('1/5 Checking input format, calls, genome build and strand...')
    genome = load_genome_data(str(input_path), assembly=args.assembly, strand=args.strand)
    snapshot = load_snapshot(args.snapshot_path) if args.snapshot_path else None
    if args.snapshot_path and snapshot is None:
        raise FileNotFoundError(f'Evidence snapshot not found: {args.snapshot_path}')
    report = build_analysis_report(genome, snapshot, markers=None if args.include_legacy_catalog else [])
    print('2/5 Checking database provenance and matching measured variants locally...')
    report['clinical'] = _get_clinical(args, genome)
    if getattr(args, 'synthetic_demo', False):
        for finding in report['clinical'].get('findings', []):
            finding['source_url'] = None
            finding['synthetic_fixture'] = True
    print('3/5 Adding ClinGen gene-disease context...')
    report['gene_validity'] = _get_gene_validity(args, report['clinical'])
    print('4/5 Checking what can be estimated from these data...')
    report['prediction_status'] = {
        'personal_disease_probabilities': 'not_estimated',
        'disease_order_over_time': 'not_estimated',
        'gene_interaction_machine_learning': 'not_trained',
        'polygenic_scores': 'not_computed',
        'reason': 'A single raw genotype file has no labeled outcomes or follow-up cohort. No externally validated, ancestry-appropriate, calibrated model has been supplied.',
        'existing_pipeline': 'https://pgsc-calc.readthedocs.io/en/latest/',
    }
    metadata = report.setdefault('metadata', {})
    metadata.update({'software_version': VERSION, 'source_code_sha256': _code_hash(),
                     'input_sha256': file_sha256(input_path), 'input_name': input_path.name,
                     'input_format': asdict(genome.stats), 'assembly': genome.assembly,
                     'assembly_source': genome.assembly_source, 'strand': genome.strand,
                     'strand_source': genome.strand_source, 'clinical_validation': 'not_clinically_validated',
                     'educational_use_only': True, 'python_version': platform.python_version(),
                     'legacy_catalog_included': args.include_legacy_catalog,
                     'snapshot_file_sha256': file_sha256(Path(args.snapshot_path)) if args.snapshot_path else None,
                     'synthetic_demo': getattr(args, 'synthetic_demo', False)})
    metadata['pipeline_steps'] = [
        {'step': 'input_validation', 'status': 'completed'},
        {'step': 'clinvar_matching', 'status': report['clinical'].get('status')},
        {'step': 'clingen_context', 'status': report['gene_validity'].get('status')},
        {'step': 'risk_and_model_readiness', 'status': 'completed_predictions_withheld'},
        {'step': 'report_and_provenance', 'status': 'completed'},
    ]
    print('5/5 Writing readable results and reproducibility manifest...')
    json_path, md_path, manifest_path = outputs
    write_json_report(report, str(json_path))
    write_markdown_report(report, str(md_path))
    manifest = {'schema_version': 1, 'metadata': metadata,
                'clinical_provenance': report['clinical'].get('database', report['clinical'].get('provenance', {})),
                'gene_validity_provenance': report['gene_validity'],
                'outputs': {p.name: file_sha256(p) for p in [json_path, md_path]}}
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    print_console_summary(report, top_n_traits=args.top_traits)
    for label, path in [('Readable report', md_path), ('Full results', json_path), ('Run manifest', manifest_path)]:
        print(f'{label}: {path.resolve()}')
    return 0


def cmd_update_evidence(args) -> int:
    from .evidence.clinical import refresh_clinvar_database
    from .evidence.gene_validity import refresh_gene_validity
    manifest = refresh_clinvar_database(assembly=args.assembly, cache_dir=args.cache_dir, source_path=args.source_vcf, show_progress=True)
    genes = refresh_gene_validity(cache_dir=args.cache_dir, source_path=args.clingen_csv)
    print(json.dumps({'clinvar': manifest, 'clingen': genes.get('metadata', {})}, indent=2))
    return 0


def cmd_status(args) -> int:
    from .evidence.clinical import load_clinvar_manifest
    from .evidence.gene_validity import load_gene_validity
    clinical = load_clinvar_manifest(assembly=args.assembly, cache_dir=args.cache_dir)
    genes = load_gene_validity(cache_dir=args.cache_dir)
    print(json.dumps({'clinvar': clinical, 'clingen': genes.get('metadata') if genes else None}, indent=2))
    return 0


def cmd_demo(args) -> int:
    from .demo import run_demo
    return run_demo(args.output_dir)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog='gene-analysis', description='Local genotype evidence annotation. No invented disease probabilities.')
    parser.add_argument('--version', action='version', version=VERSION)
    commands = parser.add_subparsers(dest='command', required=True)
    analyze = commands.add_parser('analyze', help='Run input checks, ClinVar matching, gene context and reports.')
    analyze.add_argument('genome_file', help='Raw 4/5-column genotype text or .txt.gz.')
    analyze.add_argument('--assembly', choices=['GRCh37', 'GRCh38'], help='Confirm provider genome build if absent from the header.')
    analyze.add_argument('--strand', choices=['plus', 'minus'], help='Confirm provider strand if absent from the header; clinical matching requires plus.')
    analyze.add_argument('--auto-update', action=argparse.BooleanOptionalAction, default=True)
    analyze.add_argument('--offline', action='store_true', help='Never access the network; use cached or pinned evidence.')
    analyze.add_argument('--cache-dir', help='Directory containing public databases.')
    analyze.add_argument('--clinvar-manifest', help='Use this exact ClinVar manifest; never refresh it.')
    analyze.add_argument('--clingen-snapshot', help='Use this exact ClinGen JSON snapshot; never refresh it.')
    analyze.add_argument('--snapshot-path', help='Optional historical GWAS/PGS evidence JSON; never refreshed.')
    analyze.add_argument('--include-legacy-catalog', action='store_true', help='Include the old, unverified handpicked marker panel as supplementary evidence.')
    analyze.add_argument('--stale-after-days', type=int, default=7)
    analyze.add_argument('--output-dir', default='analysis_reports')
    analyze.add_argument('--output-prefix')
    analyze.add_argument('--json-output')
    analyze.add_argument('--markdown-output')
    analyze.add_argument('--top-traits', type=int, default=10)
    analyze.set_defaults(func=cmd_analyze)
    update = commands.add_parser('update-evidence', help='Refresh official ClinVar and ClinGen public data.')
    update.add_argument('--assembly', choices=['GRCh37', 'GRCh38'], default='GRCh37')
    update.add_argument('--cache-dir')
    update.add_argument('--source-vcf', help='Import a local ClinVar VCF instead of downloading it.')
    update.add_argument('--clingen-csv', help='Import a local ClinGen CSV instead of downloading it.')
    update.set_defaults(func=cmd_update_evidence)
    status = commands.add_parser('database-status', help='Show verified local database versions and provenance.')
    status.add_argument('--assembly', choices=['GRCh37', 'GRCh38'], default='GRCh37')
    status.add_argument('--cache-dir')
    status.set_defaults(func=cmd_status)
    demo = commands.add_parser('demo', help='Run the whole annotation flow offline with tiny synthetic fixtures.')
    demo.add_argument('--output-dir', default='analysis_reports/demo')
    demo.set_defaults(func=cmd_demo)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, 'stale_after_days', 1) < 1:
        parser.error('--stale-after-days must be positive.')
    try:
        return int(args.func(args))
    except (OSError, ValueError, RuntimeError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
