"""Conservative evidence matching; this module does not estimate personal risk.

An association effect, AUROC, or R-squared cannot be converted into a person's
probability, prediction error, or heritability. Deprecated numerical helpers
therefore return ``None``. They remain importable for older callers.
"""

from __future__ import annotations

import math
import re
from typing import Any, Dict, List, Optional, Tuple


def _safe_float(value: Any) -> Optional[float]:
    """Accept an actual finite number, never extract one from an arbitrary label."""
    if isinstance(value, bool) or value is None:
        return None
    try:
        parsed = float(value)
    except (ValueError, TypeError, OverflowError):
        return None
    return parsed if math.isfinite(parsed) else None


def normalized_genotype(genotype: Optional[str]) -> Optional[str]:
    """Return a literal one- or two-base call; do not infer alleles or ploidy."""
    if not isinstance(genotype, str):
        return None
    text = genotype.strip().upper()
    if not re.fullmatch(r"[ACGT](?:[/|]?[ACGT])?", text):
        return None
    return text.replace("/", "").replace("|", "")


def allele_count(genotype: Optional[str], risk_allele: Optional[str]) -> Optional[int]:
    """Count literal bases, not dosage after strand/build harmonization.

    Missing/unsupported calls and non-single-base alleles are unknown, not zero.
    A zero is returned only for an observed supported call lacking that base.
    """
    call = normalized_genotype(genotype)
    if call is None or not isinstance(risk_allele, str):
        return None
    allele = risk_allele.strip().upper()
    if allele not in {"A", "C", "G", "T"}:
        return None
    return call.count(allele)


def _normalize_trait(value: str) -> str:
    return " ".join(re.sub(r"[^a-z0-9]+", " ", value.casefold()).split())


def association_matches_trait(association: Dict[str, Any], trait_name: str) -> bool:
    """Require a complete reported/EFO label match, never token overlap.

    This intentionally misses synonyms until a curated ontology mapping exists.
    Similar names such as type 1 and type 2 diabetes must not be conflated.
    """
    target = _normalize_trait(trait_name)
    if not target:
        return False
    labels: List[str] = []
    reported = association.get("reported_trait", [])
    if isinstance(reported, str):
        labels.append(reported)
    elif isinstance(reported, list):
        labels.extend(item for item in reported if isinstance(item, str))
    efo = association.get("efo_traits", [])
    if isinstance(efo, list):
        for item in efo:
            label = item.get("efo_trait") if isinstance(item, dict) else item
            if isinstance(label, str):
                labels.append(label)
    return any(_normalize_trait(label) == target for label in labels)


def association_effect_allele(association: Dict[str, Any], rsid: str) -> Optional[str]:
    """Accept only an explicitly identified, single-variant base allele.

    Haplotype, interaction, multi-SNP, unknown, and indel effects need a separate
    harmonization procedure and are not assigned to a single marker here.
    """
    alleles = association.get("snp_allele")
    if not isinstance(alleles, list) or len(alleles) != 1:
        return None
    item = alleles[0]
    if not isinstance(item, dict) or str(item.get("rs_id", "")).casefold() != rsid.casefold():
        return None
    allele = str(item.get("effect_allele", "")).strip().upper()
    return allele if allele in {"A", "C", "G", "T"} else None


def pick_best_association_for_trait(
    associations: List[Dict[str, Any]],
    trait_name: str,
    rsid: Optional[str] = None,
    effect_allele: Optional[str] = None,
) -> Optional[Dict[str, Any]]:
    """Select a representative exact-label association for display only.

    When rsid is supplied, also require an unambiguous allele and, if supplied,
    agreement with the catalog allele. No strand complement is guessed. The
    smallest reported p-value only selects a display record: it is not a model
    selection or a claim that the association is replicated or clinically valid.
    """
    matches = []
    for assoc in associations:
        if not isinstance(assoc, dict) or not association_matches_trait(assoc, trait_name):
            continue
        if rsid is not None:
            allele = association_effect_allele(assoc, rsid)
            if allele is None:
                continue
            if effect_allele and allele != effect_allele.strip().upper():
                continue
        elif effect_allele is not None:
            # An allele without its variant identity is insufficient evidence.
            continue
        matches.append(assoc)

    def order(item: Dict[str, Any]) -> Tuple[float, str]:
        p = _safe_float(item.get("p_value"))
        return (p if p is not None and 0 <= p <= 1 else math.inf, str(item.get("association_id", "")))

    return min(matches, key=order) if matches else None


def relative_probability_from_effect(
    allele_count_value: Optional[int],
    effect_or: Optional[float] = None,
    beta: Optional[float] = None,
    weight: Optional[float] = None,
) -> None:
    """Deprecated: association effects do not identify personal probabilities."""
    return None


def snp_error_rate(
    ci_lower: Optional[float],
    ci_upper: Optional[float],
    effect_or: Optional[float],
    beta: Optional[float],
) -> None:
    """Deprecated: a study effect CI is not an individual prediction error rate."""
    return None


def ancestry_bias_penalty(ancestry_distribution: Dict[str, Any]) -> Tuple[None, str]:
    """Describe availability without assigning an unvalidated numeric penalty."""
    if ancestry_distribution:
        note = "Study ancestry metadata is available; transferability to this person is not established."
    else:
        note = "Study ancestry coverage is unavailable; transferability is unknown."
    return None, note


def confidence_from_error_and_bias(
    error_rate: Optional[float], bias_penalty: Optional[float], coverage_penalty: float = 0.0
) -> None:
    """Deprecated: confidence requires a defined, validated statistical target."""
    return None


def pgs_trait_model_metrics(trait_model: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Preserve published model metadata without pretending to calculate a PGS."""
    available = bool(trait_model and not trait_model.get("error") and not trait_model.get("skipped"))
    performance = (trait_model or {}).get("performance") or {}
    return {
        "status": "metadata_only" if available else "unavailable",
        "personal_score": None,
        "prs_probability": None,
        "prs_error_rate": None,
        "heritability_proxy": None,
        "prs_metric": None,
        "published_performance": performance if isinstance(performance, dict) else {},
        "note": (
            "PGS Catalog records describe published models. No scoring weights were applied to this genome; "
            "published AUROC, odds ratios, and R-squared are not personal probabilities, error rates, or heritability."
        ),
    }
