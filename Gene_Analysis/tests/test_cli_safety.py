"""Network isolation, replay and user-visible failure behavior."""
import json
from pathlib import Path
import urllib.request

from Gene_Analysis.cli import main


def _forbid_network(*args, **kwargs):
    raise AssertionError('Offline operation attempted a network request')


def test_complete_offline_demo_and_repeatable_evidence(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, 'urlopen', _forbid_network)
    for name in ['first', 'second']:
        assert main(['demo', '--output-dir', str(tmp_path / name)]) == 0
    reports = []
    for name in ['first', 'second']:
        report_path = next(p for p in (tmp_path / name).glob('gene_report_*.json') if '.manifest.' not in p.name)
        report = json.loads(report_path.read_text())
        reports.append(report)
        assert report['metadata']['synthetic_demo'] is True
        assert report['clinical']['status'] == 'complete'
        assert report['clinical']['summary']['matched_variants'] == 2
        assert report['prediction_status']['gene_interaction_machine_learning'] == 'not_trained'
        assert report['clinical']['findings'][0]['absolute_risk'] is None
        assert report['clinical']['findings'][0]['source_url'] is None
        assert report['rows'] == []  # no unverified handpicked panel by default
        assert 'SYNTHETIC' in report_path.with_suffix('.md').read_text()
        manifest = json.loads(report_path.with_suffix('.manifest.json').read_text())
        assert len(manifest['metadata']['input_sha256']) == 64
        assert manifest['outputs'][report_path.name]
    assert reports[0]['clinical']['findings'] == reports[1]['clinical']['findings']
    assert reports[0]['metadata']['input_sha256'] == reports[1]['metadata']['input_sha256']


def test_unknown_build_is_withheld_not_negative(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, 'urlopen', _forbid_network)
    genome = tmp_path / 'input.txt'
    genome.write_text('rs1\t1\t100\tAG\n')
    assert main(['analyze', str(genome), '--offline', '--cache-dir', str(tmp_path / 'empty'),
                 '--output-dir', str(tmp_path), '--output-prefix', 'unknown']) == 0
    report = json.loads((tmp_path / 'unknown.json').read_text())
    assert report['clinical']['status'] == 'not_analyzed'
    assert 'unknown' in report['clinical']['reason'].lower()


def test_pinned_missing_snapshot_is_not_replaced(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, 'urlopen', _forbid_network)
    genome = tmp_path / 'input.txt'
    genome.write_text('# genome build GRCh37\n# strand plus\nrs1\t1\t100\tAG\n')
    assert main(['analyze', str(genome), '--clinvar-manifest', str(tmp_path / 'missing.json'),
                 '--output-dir', str(tmp_path)]) == 2
    assert not list(tmp_path.glob('gene_report*'))


def test_no_auto_update_without_cache_does_not_download(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, 'urlopen', _forbid_network)
    genome = tmp_path / 'input.txt'
    genome.write_text('# genome build GRCh37\n# strand plus\nrs1\t1\t100\tAG\n')
    assert main(['analyze', str(genome), '--no-auto-update', '--cache-dir', str(tmp_path / 'empty'),
                 '--output-dir', str(tmp_path)]) == 2


def test_output_cannot_overwrite_input(tmp_path):
    genome = tmp_path / 'input.txt'
    original = 'rs1\t1\t100\tAG\n'
    genome.write_text(original)
    assert main(['analyze', str(genome), '--offline', '--json-output', str(genome)]) == 2
    assert genome.read_text() == original


def test_malformed_pinned_manifest_fails_even_when_input_build_unknown(tmp_path, monkeypatch):
    monkeypatch.setattr(urllib.request, 'urlopen', _forbid_network)
    genome = tmp_path / 'input.txt'
    genome.write_text('rs1\t1\t100\tAG\n')
    manifest = tmp_path / 'bad.json'
    manifest.write_text('{"status":"incomplete"}')
    assert main(['analyze', str(genome), '--offline', '--clinvar-manifest', str(manifest),
                 '--output-dir', str(tmp_path)]) == 2
    assert not list(tmp_path.glob('gene_report*'))
