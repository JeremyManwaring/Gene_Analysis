"""Build the combined .txt (most-accurate merge of all report runs) and a
disease-likelihood dataset used by the cluster visualization."""

from __future__ import annotations

if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")



import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List

from synapse_crosscheck import cross_check, main as crosscheck_main, REPO_ROOT, load_json, SNAPSHOT_PATH
from realworld import real_world_probability

OUT_TXT = REPO_ROOT / "analysis_reports" / "combined_crosschecked_results.txt"
OUT_JSON = REPO_ROOT / "analysis_reports" / "disease_likelihood.json"


def _load_allele_counts():
    """trait -> (rsid, allele_count) from the report rows (disease markers)."""
    snap = load_json(SNAPSHOT_PATH) or {}
    markers = snap.get("markers", {})
    amap = {}
    for p in sorted((REPO_ROOT / "analysis_reports").glob("gene_report_*.json")):
        rep = load_json(p)
        if not rep:
            continue
        for r in rep.get("rows", []):
            if r.get("category") == "disease" and r.get("trait") not in amap:
                amap[r["trait"]] = (r.get("rsid"), int(r.get("risk_allele_count", 0) or 0))
    return amap, markers


def fmt(v, suffix="%"):
    return f"{v:.2f}{suffix}" if isinstance(v, (int, float)) else "n/a"


def build() -> None:
    res = crosscheck_main()
    verdicts = res["verdicts"]
    amap, markers = _load_allele_counts()

    # Compute the real-world Probability (gated cross-check) for disease traits.
    rw_by_trait = {}
    for v in verdicts:
        if v.category != "disease":
            continue
        rsid, k = amap.get(v.trait, (None, 0))
        rw_by_trait[v.trait] = real_world_probability(v.trait, k, markers.get(rsid, {}))

    # ----- combined .txt -----
    lines: List[str] = []
    lines.append("=" * 92)
    lines.append("COMBINED CROSS-CHECKED GENE ANALYSIS")
    lines.append("=" * 92)
    lines.append(f"Generated:            {datetime.now().isoformat(timespec='seconds')}")
    lines.append(f"Report runs merged:   {', '.join(res['runs'])}")
    lines.append(f"Traits combined:      {res['n_traits']}")
    lines.append(f"Evidence markers:     {res['snapshot_markers']} (cached GWAS/dbSNP/ClinVar/gnomAD)")
    lines.append(f"Synapse connector:    {res['synapse_status']}")
    lines.append("")
    lines.append("COLUMN MEANINGS")
    lines.append("  Genetic lean : raw direction/strength of your genotype for the trait (0-100).")
    lines.append("                 NOT a real-world chance; 50%/5% are neutral/absent defaults.")
    lines.append("  Likelihood   : LOW-CONFIDENCE estimate = Genetic lean weighted by evidence")
    lines.append("                 confidence. A rough comparative score, not real-world risk.")
    lines.append("  Probability  : real-world absolute-risk estimate. Filled in ONLY when the")
    lines.append("                 marker is cross-listed in >=2 databanks AND a trait-matched")
    lines.append("                 effect size exists AND a baseline prevalence is known.")
    lines.append("                 = baseline prevalence x your genotype's multiplier.")
    lines.append("                 Blank = not enough cross-checked evidence to estimate.")
    lines.append("")
    lines.append("Merge rule: per trait, the run with the HIGHEST confidence (tie: lowest error) kept.")
    lines.append("EDUCATIONAL / RESEARCH USE ONLY - not a medical diagnosis or clinical risk.")
    lines.append("=" * 92)
    lines.append("")

    # group by category
    by_cat: Dict[str, list] = {}
    for v in verdicts:
        by_cat.setdefault(v.category, []).append(v)

    for cat in sorted(by_cat):
        rows = by_cat[cat]
        lines.append(f"\n### {cat.upper()}  ({len(rows)} traits)")
        lines.append("-" * 92)
        lines.append(f"{'Trait':<48}{'GeneticLean':>12}{'Likelihood':>12}{'Probability':>13}")
        lines.append(f"{'':<48}{'(raw)':>12}{'(low-conf)':>12}{'(real-world)':>13}")
        lines.append("-" * 92)
        for v in rows:
            rw = rw_by_trait.get(v.trait, {})
            prob = rw.get("probability")
            lines.append(
                f"{v.trait[:47]:<48}"
                f"{fmt(v.relative_probability):>12}"
                f"{fmt(likelihood_score(v) * 100):>12}"
                f"{(fmt(prob) if prob is not None else '—'):>13}"
            )
            meta = f"     src: {v.chosen_run}  |  {v.note}"
            if v.evidence_sources:
                meta += f"  |  evidence: {', '.join(v.evidence_sources)}"
            lines.append(meta)
            if v.category == "disease" and prob is None and rw.get("reason"):
                lines.append(f"     probability blank: {rw['reason']}")
        lines.append("")

    # ---- Disease ranking: sort by real-world Probability first, then likelihood ----
    diseases = [v for v in verdicts if v.category == "disease"]
    for v in diseases:
        v.__dict__["likelihood"] = likelihood_score(v)
        v.__dict__["rw"] = rw_by_trait.get(v.trait, {})

    def sort_key(v):
        p = v.__dict__["rw"].get("probability")
        return (1 if p is not None else 0, p if p is not None else v.__dict__["likelihood"] * 100)

    diseases.sort(key=sort_key, reverse=True)

    lines.append("=" * 92)
    lines.append("DISEASE RANKING  (real-world Probability where available, else low-conf Likelihood)")
    lines.append("=" * 92)
    lines.append(f"{'Disease':<46}{'Probability':>13}{'Likelihood':>12}{'GeneticLean':>12}  Basis")
    lines.append("-" * 92)
    for v in diseases:
        rw = v.__dict__["rw"]
        prob = rw.get("probability")
        basis = "cross-checked real-world" if prob is not None else "low-confidence only"
        lines.append(
            f"{v.trait[:45]:<46}"
            f"{(fmt(prob) if prob is not None else '—'):>13}"
            f"{fmt(v.__dict__['likelihood'] * 100):>12}"
            f"{fmt(v.relative_probability):>12}  {basis}"
        )
    lines.append("")
    lines.append("How Probability was computed for the cross-checked rows:")
    lines.append("  risk = baseline_prevalence x (OR^your_alleles / population_average_multiplier)")
    lines.append("  Because the average person also carries risk alleles, not carrying them")
    lines.append("  usually puts you BELOW the population baseline - so many numbers went down.")
    lines.append("  These are approximate, single-marker estimates, NOT validated clinical risk.")

    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")

    # ----- JSON for the visualization -----
    payload = []
    for v in diseases:
        rw = v.__dict__["rw"]
        payload.append({
            "trait": v.trait,
            "group": disease_group(v.trait),
            "genetic_lean": v.relative_probability,
            "likelihood_low_conf": round(v.__dict__["likelihood"], 4),
            "probability_real_world": rw.get("probability"),
            "has_real_world": rw.get("probability") is not None,
            "confidence": v.confidence,
            "error_rate": v.error_rate,
            "baseline_prevalence": rw.get("baseline"),
            "or_per_allele": rw.get("or_per_allele"),
            "risk_allele_freq": rw.get("risk_allele_freq"),
            "genotype_multiplier": rw.get("multiplier"),
            "databanks": rw.get("databanks", []),
            "probability_note": rw.get("reason", ""),
            "evidence_sources": v.evidence_sources,
            "chosen_run": v.chosen_run,
        })
    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    n_rw = sum(1 for p in payload if p["has_real_world"])
    print(f"Wrote {OUT_TXT}")
    print(f"Wrote {OUT_JSON}  ({len(payload)} diseases, {n_rw} with real-world Probability)")


def likelihood_score(v) -> float:
    """Blend relative probability with confidence into a 0-1 likelihood signal."""
    p = (v.relative_probability or 0) / 100.0
    c = (v.confidence or 0) / 100.0
    # geometric-ish blend: probability weighted by how much we trust it,
    # with a floor so zero-confidence items are dimmed but not erased.
    return max(0.0, min(1.0, p * (0.35 + 0.65 * c)))


# Similarity grouping by body-system / mechanism keywords.
GROUPS = [
    ("Cancers", ["cancer", "tumor", "carcinoma", "tp53", "brca", "melanoma", "leukemia"]),
    ("Cardiometabolic", ["coronary", "artery", "cardio", "obesity", "diabetes", "cholesterol",
                          "lipid", "hypertension", "homocysteine", "thrombo", "factor v", "gout"]),
    ("Neuro / Psychiatric", ["alzheimer", "parkinson", "depression", "neurotic", "migraine",
                             "cognitive", "dementia", "schizo"]),
    ("Autoimmune / Inflammatory", ["autoimmune", "lupus", "sle", "celiac", "crohn", "multiple sclerosis",
                                    "asthma", "arthritis", "psoria", "thyroid"]),
    ("Eye / Sensory", ["macular", "degeneration", "htra1", "glaucoma", "vision", "hearing"]),
    ("Metabolic / Digestive", ["lactose", "hemochromatosis", "alcohol", "aldh2", "adh1b",
                               "gastric", "liver", "metabolism"]),
    ("Infectious / Immune susceptibility", ["covid", "ace2", "infection", "hepatitis"]),
]


def disease_group(trait: str) -> str:
    t = trait.lower()
    for name, kws in GROUPS:
        if any(k in t for k in kws):
            return name
    return "Other"


if __name__ == "__main__":
    build()
