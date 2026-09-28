"""Real-world absolute-risk model for disease markers.

Turns the raw "genetic lean" into an estimated real-world probability, but ONLY
when the evidence supports it:

  Cross-check gate (all three required, else Probability = None):
    1. the SNP is cross-listed in >= 2 databanks (GWAS/dbSNP/ClinVar/gnomAD),
    2. a GWAS association whose reported/EFO trait actually matches the disease
       supplies a usable effect size (odds ratio, or exp(beta) as a proxy),
    3. a curated literature baseline prevalence exists for that disease.

Absolute-risk model (per-allele multiplicative odds, HWE population reference):
    OR_a  = per-allele odds ratio  (or_value, else exp(beta))
    p     = risk-allele frequency  (gnomAD AF proxy)
    E     = (1 - p + p*OR_a)^2      # population-average multiplier
    M     = OR_a^k / E              # this genotype's multiplier vs population
    risk  = clamp(baseline * M, 0..0.95)

Because the average genotype carries risk alleles too, a person with k=0 risk
copies typically lands BELOW the baseline prevalence -- so many numbers go DOWN.
That is intended: it is the honest, evidence-anchored estimate.

Everything is approximate and for education/research only -- not clinical risk.
"""

from __future__ import annotations

if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")



import math
import re
from typing import Any, Dict, List, Optional, Tuple


def _tokenize(value: str) -> List[str]:
    return [t for t in re.sub(r"[^a-z0-9]+", " ", value.lower()).split() if t]


def _safe_float(value: Any) -> Optional[float]:
    if isinstance(value, (int, float)):
        return float(value)
    return None


def pick_best_association_for_trait(
    associations: List[Dict[str, Any]], trait_name: str
) -> Optional[Dict[str, Any]]:
    """Rank associations by token overlap with the disease name, then by p-value."""
    if not associations:
        return None
    trait_tokens = set(_tokenize(trait_name))
    ranked: List[Tuple[float, float, Dict[str, Any]]] = []
    for assoc in associations:
        chunks: List[str] = []
        rt = assoc.get("reported_trait", [])
        chunks += [str(x) for x in rt] if isinstance(rt, list) else [str(rt)]
        for item in assoc.get("efo_traits", []) or []:
            chunks.append(str(item.get("efo_trait", "")) if isinstance(item, dict) else str(item))
        assoc_tokens = set(_tokenize(" ".join(chunks)))
        overlap = len(trait_tokens.intersection(assoc_tokens))
        p_val = _safe_float(assoc.get("p_value"))
        p_rank = p_val if p_val is not None else 10.0
        ranked.append((float(overlap), -p_rank, assoc))
    ranked.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return ranked[0][2]


# Curated literature baseline lifetime/period prevalences (general adult pop.,
# predominantly European-ancestry reference where the markers were derived).
# Values are approximate order-of-magnitude figures from public epidemiology.
BASELINE_PREVALENCE: Dict[str, float] = {
    "Celiac disease": 0.010,
    "Coronary artery disease": 0.300,
    "Crohn's disease": 0.003,
    "Obesity": 0.400,
    "Parkinson's disease": 0.020,
    "Psoriasis": 0.030,
    "Prostate cancer": 0.130,
    "Depression / neuroticism": 0.200,
    "Migraine": 0.150,
    "Lupus (SLE)": 0.0010,
    "Gout (ABCG2)": 0.040,
    "Age-related macular degeneration (HTRA1)": 0.100,
    "Gastric cancer (LINC00673)": 0.008,
    "Hemochromatosis": 0.004,
    "Lactose intolerance": 0.170,
}

# Traits deliberately excluded from real-world probability because they are
# biomarkers / phenotypes / deterministic reactions, not disease outcomes,
# or too non-specific to anchor a prevalence to.
NON_DISEASE_TRAITS = {
    "Homocysteine metabolism (MTHFR)",
    "Alcohol metabolism (ADH1B)",
    "Alcohol flush reaction (ALDH2 deficiency)",
    "Skin pigmentation (OCA2)",
    "Cancer risk (TP53)",
    "Autoimmune disease (FCGR2A)",
}


def _databank_count(ev: Dict[str, Any]) -> List[str]:
    banks: List[str] = []
    if (ev.get("gwas_snp") or {}) and not (ev.get("gwas_snp") or {}).get("error"):
        banks.append("GWAS")
    if (ev.get("dbsnp") or {}) and not (ev.get("dbsnp") or {}).get("error"):
        banks.append("dbSNP")
    clin = ev.get("clinvar") or {}
    if clin and not clin.get("error") and clin.get("count"):
        banks.append("ClinVar")
    gno = ev.get("gnomad") or {}
    if gno and not gno.get("error") and not gno.get("skipped"):
        banks.append("gnomAD")
    return banks


def _risk_allele_freq(ev: Dict[str, Any]) -> Optional[float]:
    gno = ev.get("gnomad") or {}
    af = gno.get("genome_af")
    if af is None:
        af = gno.get("exome_af")
    if isinstance(af, (int, float)) and 0 < af < 1:
        return float(af)
    return None


def _matched_effect(trait: str, associations: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Return the trait-matched association only if it actually overlaps the
    disease name AND carries an effect size. Otherwise None (fails cross-check)."""
    if not associations:
        return None
    best = pick_best_association_for_trait(associations, trait)
    if not best:
        return None
    # require genuine token overlap between disease name and association trait
    trait_tokens = set(_tokenize(trait))
    chunks: List[str] = []
    rt = best.get("reported_trait", [])
    chunks += [str(x) for x in rt] if isinstance(rt, list) else [str(rt)]
    for item in best.get("efo_traits", []) or []:
        chunks.append(str(item.get("efo_trait", "")) if isinstance(item, dict) else str(item))
    assoc_tokens = set(_tokenize(" ".join(chunks)))
    if not trait_tokens.intersection(assoc_tokens):
        return None
    if best.get("or_value") or best.get("beta") is not None:
        return best
    return None


def per_allele_or(assoc: Dict[str, Any]) -> Optional[float]:
    orv = assoc.get("or_value")
    if isinstance(orv, (int, float)) and orv > 0:
        return float(orv)
    beta = assoc.get("beta")
    if isinstance(beta, (int, float)):
        return math.exp(float(beta))  # proxy: log-odds -> odds ratio
    return None


def real_world_probability(trait: str, allele_count: int, ev: Dict[str, Any]) -> Dict[str, Any]:
    """Compute the gated real-world probability estimate for one disease marker."""
    result: Dict[str, Any] = {
        "probability": None,
        "baseline": None,
        "or_per_allele": None,
        "risk_allele_freq": None,
        "multiplier": None,
        "databanks": _databank_count(ev),
        "reason": "",
    }

    if trait in NON_DISEASE_TRAITS:
        result["reason"] = "biomarker/phenotype, not a disease outcome"
        return result

    baseline = BASELINE_PREVALENCE.get(trait)
    if baseline is None:
        result["reason"] = "no curated baseline prevalence"
        return result
    result["baseline"] = baseline

    banks = result["databanks"]
    if len(banks) < 2:
        result["reason"] = f"not cross-listed (only {banks or 'none'})"
        return result

    assoc = _matched_effect(trait, ev.get("gwas_associations", []) or [])
    if not assoc:
        result["reason"] = "no trait-matched effect size in databanks"
        return result

    or_a = per_allele_or(assoc)
    if not or_a or or_a <= 0:
        result["reason"] = "effect size unusable"
        return result
    result["or_per_allele"] = round(or_a, 4)

    p = _risk_allele_freq(ev)
    if p is None:
        # fall back to association's own risk_frequency if present
        rf = assoc.get("risk_frequency")
        p = float(rf) if isinstance(rf, (int, float)) and 0 < rf < 1 else 0.2
    result["risk_allele_freq"] = round(p, 4)

    k = max(0, min(2, int(allele_count or 0)))
    pop_mult = (1 - p + p * or_a) ** 2  # expected multiplier over population
    if pop_mult <= 0:
        result["reason"] = "degenerate population frequency"
        return result
    mult = (or_a ** k) / pop_mult
    result["multiplier"] = round(mult, 3)

    risk = max(0.0, min(0.95, baseline * mult))
    result["probability"] = round(risk * 100.0, 2)  # percent
    result["reason"] = "cross-checked: baseline x genotype multiplier"
    return result
