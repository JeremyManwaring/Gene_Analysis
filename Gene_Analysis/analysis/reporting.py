"""Readable clinical-source annotations with explicit scope and provenance."""

from __future__ import annotations

import html
import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List


MAX_MARKDOWN_FINDINGS = 100


def _ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def write_json_report(report: Dict[str, Any], path: str) -> str:
    output_path = Path(path)
    _ensure_parent(output_path)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, allow_nan=False)
        handle.write("\n")
    return str(output_path)


def _cell(value: Any) -> str:
    if value is None or value == "":
        return "not available"
    if isinstance(value, (dict, list)):
        value = json.dumps(value, ensure_ascii=False, sort_keys=True)
    return html.escape(str(value)).replace("|", "\\|").replace("\r", " ").replace("\n", " ")


def _display(value: Any) -> str:
    return str(value or "not available").replace("_", " ")


def _clinical_lines(report: Dict[str, Any]) -> List[str]:
    clinical = report.get("clinical") or {}
    lines = ["## ClinVar variant annotations", ""]
    if clinical.get("status") != "complete":
        lines += [
            "**Clinical matching was not completed.** " + _cell(clinical.get("reason") or "No clinical annotation results were supplied."),
            "", "No result here should be treated as a negative clinical test.", "",
        ]
        return lines
    summary = clinical.get("summary") or {}
    lines += [
        f"Input records checked: **{summary.get('input_records', 0)}**. "
        f"Matched alternate variant records: **{summary.get('matched_variants', 0)}**. "
        f"Flagged for professional review: **{summary.get('professional_review_candidates', 0)}**.",
        "",
        "These are ClinVar's submitted variant classifications. They do not establish that you have, "
        "will develop, or carry an inherited condition. Confirm any medically relevant finding with a qualified clinician and a clinical test.",
        "",
        f"Expert panel / practice guideline records: {summary.get('expert_reviewed_matches', 0)}. "
        f"Records with conflicting classifications: {summary.get('conflicting_matches', 0)}.",
        "",
    ]
    findings = clinical.get("findings") or []
    if findings:
        # Prioritize review workflow, never calculated disease likelihood.
        ordered = sorted(findings, key=lambda item: (
            not bool(item.get("needs_professional_review")),
            not bool(item.get("conflicting_classifications")),
            -int(item.get("review_stars") or 0),
            str(item.get("chromosome", "")), int(item.get("position") or 0),
        ))
        visible = ordered[:MAX_MARKDOWN_FINDINGS]
        lines += [
            f"Showing {len(visible)} of {len(findings)} matched variant records, with records needing review first. "
            "The JSON report contains every match and its full source fields. This order is not a disease-risk ranking.",
            "",
            "| Gene | Variant / genotype | ClinVar classification | Review status | Source condition(s) | Record |",
            "|---|---|---|---|---|---|",
        ]
        for finding in visible:
            variant = (
                f"{finding.get('chromosome', '?')}:{finding.get('position', '?')} "
                f"{finding.get('reference', '?')}>{finding.get('alternate', '?')}; "
                f"call {finding.get('genotype', '?')}"
            )
            stars = finding.get("review_stars", 0)
            review = f"{stars}/4 stars; {_display(finding.get('review_status'))}"
            variation_id = str(finding.get("clinvar_variation_id", ""))
            source = (
                f"[ClinVar {variation_id}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{variation_id}/)"
                if variation_id.isdigit() else "not available"
            )
            if (report.get('metadata') or {}).get('synthetic_demo'):
                source = 'Synthetic fixture; no live record'
            lines.append("| " + " | ".join([
                _cell(finding.get("gene")), _cell(variant),
                _cell(finding.get("germline_classification")), _cell(review),
                _cell(finding.get("conditions")), source,
            ]) + " |")
        lines += [
            "", "Stars describe ClinVar review status (3: expert panel; 4: practice guideline). "
            "They are not a personal confidence score. Aggregate classifications can cover multiple conditions; "
            "open the source record to review each assertion, inheritance, and conflicts.", "",
        ]
    else:
        lines += ["No supported observed alternate alleles matched this ClinVar snapshot. This does not rule out any condition.", ""]
    lines += [
        "### Coverage limits", "",
        "| Check | Count |", "|---|---:|",
    ]
    for key, label in (
        ("records_with_match", "Input records with a supported alternate allele match"),
        ("coordinates_without_clinvar", "Input coordinates without a ClinVar record"),
        ("reference_only_records", "Input records matching only the reference base at supported sites"),
        ("unsupported_genotypes", "Unsupported calls"),
        ("invalid_coordinates", "Invalid coordinates"),
        ("unsupported_variant_records", "ClinVar candidate records outside SNV support"),
        ("allele_mismatches", "ClinVar candidate records with incompatible alleles"),
    ):
        if key in summary:
            lines.append(f"| {label} | {_cell(summary[key])} |")
    lines += [
        "", "Counts mix input records and candidate ClinVar records as labeled, so they should not be added together.",
        "", _cell(clinical.get("scope") or "Only supported measured alleles are evaluated."), "",
    ]
    return lines


def _provenance_lines(report: Dict[str, Any]) -> List[str]:
    meta = report.get("metadata") or {}
    clinical = report.get("clinical") or {}
    database = clinical.get("database") or clinical.get("provenance") or {}
    lines = ["## Reproduce this report", "", "The JSON report and run manifest retain the full evidence and provenance.", ""]
    for key, label in (
        ("software_version", "Software version"), ("python_version", "Python version"),
        ("input_sha256", "Input SHA-256"), ("source_code_sha256", "Analysis code SHA-256"),
        ("snapshot_file_sha256", "Supplementary snapshot SHA-256"),
    ):
        if meta.get(key):
            lines.append(f"- {label}: `{_cell(meta[key])}`")
    for key, label in (
        ("source_release", "ClinVar source release"), ("fetched_at", "ClinVar retrieved at"),
        ("source_kind", "ClinVar source kind"), ("source_sha256", "ClinVar source SHA-256"),
        ("database_sha256", "ClinVar database SHA-256"), ("manifest_sha256", "ClinVar manifest SHA-256"),
        ("variant_count", "ClinVar indexed records"), ("snv_count", "ClinVar indexed SNV records"),
    ):
        if database.get(key) is not None:
            lines.append(f"- {label}: `{_cell(database[key])}`")
    stale = clinical.get("cache_stale", database.get("stale"))
    if stale is not None:
        lines.append(f"- ClinVar cache exceeds the configured age threshold: **{'yes' if stale else 'no'}**.")
    if clinical.get("pinned"):
        lines.append("- This run used an explicitly pinned ClinVar snapshot.")
    lines += ["", "A recent retrieval does not guarantee that every submitted classification is current or correct.", ""]
    return lines


def write_markdown_report(report: Dict[str, Any], path: str) -> str:
    output_path = Path(path)
    _ensure_parent(output_path)
    meta = report.get("metadata") or {}
    lines = ["# Genome evidence report", ""]
    if meta.get("synthetic_demo"):
        lines += ["**SYNTHETIC DEMO — fabricated input and evidence for testing the workflow. Not a person's results or a live clinical database.**", ""]
    lines += [
        f"Generated: {_cell(meta.get('generated_at', datetime.now().isoformat()))}", "",
        "**This pipeline annotates measured variants with source evidence. It is not a clinically validated diagnostic or disease-prediction test.**", "",
    ]
    quality = report.get("input_quality") or {}
    stats = quality.get("stats") or meta.get("input_format") or {}
    if quality or stats:
        lines += [
            "## Input checks", "",
            f"Genome build: **{_cell(quality.get('assembly', meta.get('assembly', 'unknown')))}**. "
            f"Strand: **{_cell(quality.get('strand', meta.get('strand', 'unknown')))}**. "
            f"Usable unique calls: **{stats.get('parsed_rows', 0)}**.", "",
            f"Build evidence: {_cell(quality.get('assembly_source', meta.get('assembly_source')))}. "
            f"Strand evidence: {_cell(quality.get('strand_source', meta.get('strand_source')))}.", "",
        ]
        for key, label in (
            ("skipped_no_calls", "No-calls omitted"), ("skipped_unsupported_genotypes", "Unsupported calls omitted"),
            ("skipped_invalid_rows", "Invalid rows omitted"), ("duplicate_rows", "Duplicate rows"),
            ("excluded_conflicting_markers", "Conflicting marker calls excluded"),
        ):
            if stats.get(key):
                lines.append(f"- {label}: {_cell(stats[key])}")
        lines.append("")
    steps = meta.get("pipeline_steps") or report.get("pipeline_steps") or []
    if steps:
        lines += ["## Pipeline completion", "", "| Step | Status |", "|---|---|"]
        for step in steps:
            lines.append(f"| {_cell(_display(step.get('step')))} | {_cell(_display(step.get('status')))} |")
        lines.append("")
    lines += _clinical_lines(report)
    gene_validity = report.get("gene_validity") or {}
    if gene_validity:
        lines += ["## ClinGen gene–disease context", "", f"Status: **{_cell(_display(gene_validity.get('status')))}**.", ""]
        for key in ("reason", "scope"):
            if gene_validity.get(key):
                lines += [_cell(gene_validity[key]), ""]
        gene_meta = gene_validity.get("metadata") or {}
        gene_summary = gene_meta.get("summary") or {}
        if gene_summary:
            lines += [
                f"This export contains **{gene_summary.get('curations', 0)}** gene–disease curations across "
                f"**{gene_summary.get('genes', 0)}** genes. "
                f"**{gene_summary.get('strong_or_definitive_curations', 0)}** curations across "
                f"**{gene_summary.get('genes_with_strong_or_definitive_curation', 0)}** genes are classified Strong or Definitive.",
                "", "The export also retains Limited, Disputed, Refuted, and other classifications; "
                "these counts are database coverage, not personal risk estimates.", "",
            ]
            classification_counts = gene_summary.get("classification_counts") or {}
            if classification_counts:
                lines += ["| ClinGen classification | Curations in this export |", "|---|---:|"]
                for label, count in sorted(classification_counts.items()):
                    lines.append(f"| {_cell(label)} | {_cell(count)} |")
                lines.append("")
        findings = (report.get("clinical") or {}).get("findings") or []
        with_context = sum(bool(item.get("clingen_gene_validity")) for item in findings)
        exact_context = sum(any(
            context.get("disease_identifier_in_clinvar_aggregate")
            for context in item.get("clingen_gene_validity", [])
        ) for item in findings)
        if findings:
            lines += [
                f"Matched variant records with gene-level context: **{with_context} / {len(findings)}**. "
                f"Of these, **{exact_context}** have an exact disease identifier shared with the aggregate ClinVar record. "
                "Full curations, inheritance labels, classifications, and source links are retained per finding in JSON.", "",
            ]
        for key, label in (
            ("source_created_at", "Source export created"), ("fetched_at", "Retrieved"),
            ("source_sha256", "Source SHA-256"),
        ):
            if gene_meta.get(key) is not None:
                lines.append(f"- {label}: `{_cell(gene_meta[key])}`")
        snapshot_hash = gene_validity.get("snapshot_sha256") or (
            (report.get("clinical") or {}).get("clingen_provenance") or {}
        ).get("snapshot_sha256")
        if snapshot_hash:
            lines.append(f"- Snapshot SHA-256: `{_cell(snapshot_hash)}`")
        if gene_validity.get("cache_stale") is not None:
            lines.append(f"- Cache exceeds the configured age threshold: **{'yes' if gene_validity['cache_stale'] else 'no'}**.")
        lines += ["", "A gene–disease relationship does not prove that a particular variant causes disease. Condition, inheritance, and variant evidence must also match.", ""]
    prediction = report.get("prediction_status") or report.get("model_status") or {}
    lines += [
        "## What was not estimated", "",
        "- Personal disease probabilities or the order in which conditions might develop.",
        "- A machine-learning model of gene combinations.",
        "- A personal polygenic score, prediction error rate, confidence percentage, or heritability.", "",
        _cell(prediction.get("reason") or "No externally validated and calibrated individual prediction model is supplied."), "",
        "Published odds ratios, beta coefficients, AUROC, and R-squared are study or model statistics. "
        "They are not a person's probability of disease. PGS Catalog metadata alone is not a computed score.", "",
    ]
    traits = report.get("trait_summaries") or []
    rows = report.get("rows") or []
    if traits or rows:
        lines += [
            "## Optional historical marker catalog", "",
            "These are unverified supplementary catalog labels, separate from the ClinVar allele matches. "
            "They do not validate gene, trait, or disease interpretations. Similar trait names are not automatically equated.", "",
            f"Unique selected markers observed: {meta.get('unique_catalog_markers_found', 0)} / {meta.get('unique_catalog_markers', 0)}.", "",
            "| Category | Catalog trait | Observed calls | Exact trait and allele study matches |",
            "|---|---|---:|---:|",
        ]
        for item in sorted(traits, key=lambda value: (str(value.get("category", "")), str(value.get("trait", "")))):
            lines.append(f"| {_cell(item.get('category'))} | {_cell(item.get('trait'))} | "
                         f"{item.get('markers_found', 0)}/{item.get('markers_total', 0)} | {item.get('matched_study_records', 0)} |")
        lines += ["", "Full source records, raw calls, and unharmonized allele counts are in the JSON report. "
                  "No probability or clinical conclusion is calculated from this catalog.", ""]
    lines += _provenance_lines(report)
    lines += ["## Interpretation limits", ""]
    limits = (report.get("clinical") or {}).get("limitations") or report.get("limitations") or []
    for limit in limits:
        lines.append(f"- {_cell(limit)}")
    lines += ["- A missing call, unmatched variant, or no flagged finding is not evidence that a disease is absent."]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(output_path)


def print_console_summary(report: Dict[str, Any], top_n_traits: int = 20) -> None:
    """Summarize completion and coverage; never rank conditions by apparent risk."""
    meta = report.get("metadata") or {}
    if meta.get("synthetic_demo"):
        print("\nSYNTHETIC DEMO: fabricated data and evidence; not personal or clinical results.")
    print("\nGENOME EVIDENCE SUMMARY")
    clinical = report.get("clinical") or {}
    if clinical.get("status") == "complete":
        summary = clinical.get("summary") or {}
        database = clinical.get("database") or {}
        print(f"ClinVar source release: {database.get('source_release', 'unknown')}")
        print(f"Input records checked: {summary.get('input_records', 0)}")
        print(f"Matched alternate variant records: {summary.get('matched_variants', 0)}")
        print(f"Flagged for professional review: {summary.get('professional_review_candidates', 0)}")
        print(f"Records with conflicting classifications: {summary.get('conflicting_matches', 0)}")
        if clinical.get("cache_stale", database.get("stale")):
            print("ClinVar snapshot is older than the configured age threshold.")
    else:
        print("Clinical matching withheld: " + str(clinical.get("reason") or "no clinical annotation results supplied"))
    gene_validity = report.get("gene_validity") or {}
    if gene_validity:
        print(f"ClinGen gene-disease context: {_display(gene_validity.get('status'))}")
    traits = report.get("trait_summaries") or []
    if traits and top_n_traits > 0:
        print("\nOptional unverified catalog coverage (alphabetical):")
        for item in sorted(traits, key=lambda value: str(value.get("trait", "")))[:top_n_traits]:
            print(f"  {item.get('trait', '')}: {item.get('markers_found', 0)}/{item.get('markers_total', 0)} calls")
    print("Personal probabilities, disease timing, and gene-combination ML: not estimated.")
    print("These are source annotations, not a diagnosis. No match does not rule out disease.\n")
