"""Build auditable evidence summaries without unvalidated personal predictions."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, is_dataclass
from datetime import UTC, datetime, timedelta
from typing import Any, Dict, List, Optional, Tuple

from .catalog import MarkerDefinition, build_marker_catalog
from .scoring import (
    _safe_float,
    allele_count,
    ancestry_bias_penalty,
    association_effect_allele,
    association_matches_trait,
    normalized_genotype,
    pgs_trait_model_metrics,
    pick_best_association_for_trait,
)


UNESTIMATED = {
    "relative_probability": None,
    "error_rate": None,
    "confidence": None,
    "heritability_proxy": None,
}


def _parse_isoish(value: Optional[str]) -> Optional[datetime]:
    if not isinstance(value, str) or not value.strip():
        return None
    text = value.strip().replace("Z", "+00:00")
    try:
        dt = datetime.fromisoformat(text)
        return dt.replace(tzinfo=UTC) if dt.tzinfo is None else dt.astimezone(UTC)
    except ValueError:
        pass
    try:
        return datetime.strptime(text, "%Y/%m/%d %H:%M").replace(tzinfo=UTC)
    except ValueError:
        return None


def _valid_payload(payload: Any) -> bool:
    return isinstance(payload, dict) and bool(payload) and not any(
        payload.get(key) for key in ("error", "_error", "skipped")
    )


def _snp_type(marker_payload: Dict[str, Any]) -> str:
    gwas_snp = marker_payload.get("gwas_snp", {})
    if _valid_payload(gwas_snp):
        for key in ("functional_class", "most_severe_consequence"):
            if gwas_snp.get(key):
                return str(gwas_snp[key])
    dbsnp = marker_payload.get("dbsnp", {})
    if _valid_payload(dbsnp) and dbsnp.get("fxn_class"):
        return str(dbsnp["fxn_class"]).split(",", 1)[0]
    return "unknown"


def _source_tags(marker_payload: Dict[str, Any], trait_model: Optional[Dict[str, Any]]) -> List[str]:
    tags = []
    for key, label in (
        ("gwas_snp", "GWAS Catalog"), ("dbsnp", "dbSNP"),
        ("clinvar", "ClinVar lookup metadata"), ("gnomad", "gnomAD"),
    ):
        if _valid_payload(marker_payload.get(key)):
            tags.append(label)
    if marker_payload.get("gwas_associations"):
        tags.append("GWAS Catalog")
    if _valid_payload(trait_model):
        tags.append("PGS Catalog metadata")
    return sorted(set(tags))


def _evidence_freshness(marker_payload: Dict[str, Any], stale_after_days: int) -> str:
    # A source's publication/edit date does not establish when we checked it.
    retrieved = _parse_isoish(marker_payload.get("retrieved_at") or marker_payload.get("fetched_at"))
    if retrieved is None:
        return "unknown"
    age = datetime.now(UTC) - retrieved
    if age < timedelta(0):
        return "unknown"
    return "stale" if age > timedelta(days=stale_after_days) else "recently_retrieved"


def _out_of_date_flag(marker_payload: Dict[str, Any], stale_after_days: int) -> bool:
    """Compatibility flag: unknown retrieval time is not verified fresh."""
    return _evidence_freshness(marker_payload, stale_after_days) != "recently_retrieved"


def build_analysis_report(
    genome_data,
    evidence_snapshot: Optional[Dict[str, Any]],
    markers: Optional[List[MarkerDefinition]] = None,
) -> Dict[str, Any]:
    """Summarize supplementary catalog coverage and source records.

    Clinical allele-level annotations, if available, are added separately by the
    CLI. These hand-curated marker labels do not establish clinical validity.
    """
    markers = build_marker_catalog() if markers is None else markers
    evidence_snapshot = evidence_snapshot or {}
    snapshot_meta = evidence_snapshot.get("metadata") or {}
    try:
        stale_after_days = max(1, int(snapshot_meta.get("stale_after_days", 30)))
    except (ValueError, TypeError):
        stale_after_days = 30
    marker_evidence = evidence_snapshot.get("markers") or {}
    trait_models = evidence_snapshot.get("trait_models") or {}
    rows: List[Dict[str, Any]] = []
    grouped_rows: Dict[Tuple[str, str], List[Dict[str, Any]]] = defaultdict(list)

    for marker in markers:
        genotype = genome_data.get_genotype(marker.rsid)
        call = normalized_genotype(genotype)
        rsid_payload = marker_evidence.get(marker.rsid) or {}
        associations = rsid_payload.get("gwas_associations") or []
        if not isinstance(associations, list):
            associations = []
        exact_matches = [
            assoc for assoc in associations
            if isinstance(assoc, dict) and association_matches_trait(assoc, marker.trait)
        ]
        best_assoc = pick_best_association_for_trait(
            exact_matches, marker.trait, marker.rsid, marker.risk_allele
        )
        effect_allele = association_effect_allele(best_assoc, marker.rsid) if best_assoc else None
        if not exact_matches:
            association_status = "no_exact_trait_match"
        elif best_assoc is None:
            association_status = "effect_allele_unverified_or_mismatched"
        else:
            association_status = "study_record_matched_not_personally_interpreted"

        trait_model = trait_models.get(marker.trait)
        ancestry_dist = (((trait_model or {}).get("score") or {}).get("ancestry_distribution") or {})
        _, bias_note = ancestry_bias_penalty(ancestry_dist)
        freshness = _evidence_freshness(rsid_payload, stale_after_days)
        assoc = best_assoc or {}
        row = {
            "category": marker.category,
            "trait": marker.trait,
            "rsid": marker.rsid,
            "snp_type": _snp_type(rsid_payload),
            "genotype": genotype,
            "genotype_status": "missing" if not genotype else "observed" if call else "unsupported_call",
            "catalog_allele": marker.risk_allele,
            # A study effect allele can be protective or concern a non-disease
            # measurement; the old risk_allele label must not assert risk.
            "risk_allele": None,
            "effect_allele": effect_allele,
            "risk_allele_count": None,
            "raw_effect_allele_count": allele_count(genotype, effect_allele),
            "allele_harmonization": "not_performed",
            **UNESTIMATED,
            "prediction_status": "not_estimated",
            "association_status": association_status,
            "exact_trait_association_count": len(exact_matches),
            "association_id": assoc.get("association_id"),
            "association": best_assoc,
            "study_odds_ratio": _safe_float(assoc.get("or_value")),
            "study_beta": _safe_float(assoc.get("beta")),
            "study_p_value": _safe_float(assoc.get("p_value")),
            # The former untyped effect_size field conflated ORs and betas.
            "effect_size": None,
            "ci_lower": _safe_float(assoc.get("ci_lower")),
            "ci_upper": _safe_float(assoc.get("ci_upper")),
            "bias_note": bias_note,
            "source": ", ".join(_source_tags(rsid_payload, trait_model)),
            "evidence_freshness": freshness,
            "out_of_date": freshness != "recently_retrieved",
            "description": marker.description,
            "gene": marker.gene,
            "catalog_annotation_status": "unverified_supplementary_label",
            "clinical_interpretation": "not_performed",
        }
        rows.append(row)
        grouped_rows[(marker.category, marker.trait)].append(row)

    trait_summaries = []
    for (category, trait), trait_rows in sorted(grouped_rows.items()):
        trait_summaries.append({
            "category": category,
            "trait": trait,
            "markers_total": len(trait_rows),
            "markers_found": sum(row["genotype_status"] == "observed" for row in trait_rows),
            "matched_study_records": sum(row["association"] is not None for row in trait_rows),
            **UNESTIMATED,
            "prediction_status": "not_estimated",
            "pgs": pgs_trait_model_metrics(trait_models.get(trait)),
            "bias_note": trait_rows[0]["bias_note"],
            "sources": sorted({
                src.strip() for row in trait_rows for src in row["source"].split(",") if src.strip()
            }),
        })

    category_summaries = []
    for category in sorted({row["category"] for row in rows}):
        subset = [row for row in rows if row["category"] == category]
        category_summaries.append({
            "category": category,
            "traits": len({row["trait"] for row in subset}),
            "markers_total": len(subset),
            "markers_found": sum(row["genotype_status"] == "observed" for row in subset),
            **UNESTIMATED,
            "prediction_status": "not_estimated",
        })

    unique_markers = {row["rsid"] for row in rows}
    observed = {row["rsid"] for row in rows if row["genotype_status"] == "observed"}
    stats = getattr(genome_data, "stats", None)
    return {
        "metadata": {
            "schema_version": "2.0",
            "generated_at": datetime.now(UTC).isoformat(),
            "report_type": "evidence_annotation_not_clinical_prediction",
            "total_markers": len(rows),
            "markers_found": sum(row["genotype_status"] == "observed" for row in rows),
            "unique_catalog_markers": len(unique_markers),
            "unique_catalog_markers_found": len(observed),
            "snapshot_id": snapshot_meta.get("snapshot_id"),
            "snapshot_created_at": snapshot_meta.get("created_at"),
            "snapshot_stale_after_days": stale_after_days,
            "clinical_validation": "not_established_for_this_pipeline",
        },
        "input_quality": {
            key: getattr(genome_data, key, "unknown")
            for key in ("assembly", "assembly_source", "strand", "strand_source", "provider")
        } | {"stats": asdict(stats) if is_dataclass(stats) else {}},
        "model_status": {
            "absolute_disease_probabilities": "not_estimated",
            "longitudinal_disease_probabilities": "not_estimated",
            "gene_combination_machine_learning": "not_trained_or_validated",
            "polygenic_scores": "not_computed",
            "cross_disease_ranking": "not_supported",
            "reason": (
                "No trained, externally validated and calibrated individual prediction model is supplied. "
                "Longitudinal estimates require a defined time horizon and outcome-labeled follow-up cohort, "
                "clinical covariates, appropriate ancestry validation, and baseline incidence. "
                "A single genome and database associations cannot establish those probabilities."
            ),
        },
        "limitations": [
            "The supplementary catalog is a selected marker list, not all genes or all clinically relevant variants.",
            "A missing call or no database match does not mean absence of disease risk or a benign result.",
            "Database assertions and GWAS associations are source evidence, not a diagnosis or validated personal prediction.",
            "Supplementary catalog descriptions and alleles are historical annotations, not independently validated claims.",
            "Raw allele counts are literal input matches; study-effect interpretation requires build and strand harmonization.",
            "Model performance metadata is not a computed polygenic score. Cross-trait likelihood ranking is not supported.",
        ],
        "rows": rows,
        "trait_summaries": trait_summaries,
        "category_summaries": category_summaries,
        "source_evidence": {
            "markers": {rsid: marker_evidence[rsid] for rsid in sorted(unique_markers) if rsid in marker_evidence},
            "trait_models": {trait: trait_models[trait] for _, trait in grouped_rows if trait in trait_models},
        },
    }
