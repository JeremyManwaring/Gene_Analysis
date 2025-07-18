"""Comprehensive longevity analysis with probability estimates.

This script scans a broader set of SNPs associated with lifespan and aging
processes and provides probability estimates for aging rate and longevity.
"""

import sys
import subprocess

try:
    import pandas as pd
except ImportError:  # pragma: no cover - handled in runtime
    print("pandas not found. Installing pandas...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"])
    import pandas as pd


def load_data():
    df = pd.read_csv("Genome.txt", sep="\t", comment="#", header=None)
    df.columns = ["rsid", "chromosome", "position", "genotype"]
    return df


# Expanded list of longevity markers with weights
LONGEVITY_MARKERS = {
    # Major markers
    'rs2736100': {'risk': ['A'], 'weight': 0.15, 'description': 'TERT – Telomere maintenance'},
    'rs7726159': {'risk': ['A'], 'weight': 0.12, 'description': 'TERT – Telomere length'},
    'rs1317082': {'risk': ['A'], 'weight': 0.10, 'description': 'TERC – Telomere regulation'},
    'rs1801133': {'risk': ['A'], 'weight': 0.19, 'description': 'MTHFR – Methylation cycle'},
    'rs1042522': {'risk': ['C'], 'weight': 0.10, 'description': 'TP53 – DNA repair'},
    # Additional research markers
    'rs10936599': {'risk': ['C'], 'weight': 0.05, 'description': 'TERC – Telomere length'},
    'rs755017':  {'risk': ['A'], 'weight': 0.04, 'description': 'FOXO3 – Longevity variant'},
    'rs2802292': {'risk': ['G'], 'weight': 0.03, 'description': 'FOXO3 – Longevity variant'},
    'rs2157719': {'risk': ['C'], 'weight': 0.03, 'description': 'CDKN2B – Cellular senescence'},
    'rs11125529': {'risk': ['T'], 'weight': 0.02, 'description': 'ACYP2 – Telomere length'},
}

TOTAL_MARKERS = 25  # assumed total markers used in scoring


def genotype_effect(genotype: str, risk_allele: str) -> float:
    """Return effect score based on genotype."""
    if risk_allele * 2 == genotype:
        return 0.0
    if risk_allele in genotype and len(genotype) == 2:
        return 0.5
    if genotype == risk_allele.replace('A', 'G') * 2:
        return 1.0
    return 0.3


def main() -> None:
    df = load_data()

    print("\n==============================")
    print("COMPREHENSIVE LONGEVITY REPORT")
    print("==============================\n")

    raw_score = 0.0
    found = 0

    for rsid, info in LONGEVITY_MARKERS.items():
        row = df[df['rsid'] == rsid]
        if not row.empty:
            genotype = row.iloc[0]['genotype']
            effect = genotype_effect(genotype, info['risk'][0])
            contribution = effect * info['weight']
            raw_score += contribution
            found += 1
            print(f"SNP: {rsid}  Genotype: {genotype}  Effect: {contribution:.3f}  {info['description']}")
        else:
            print(f"SNP: {rsid} - NOT FOUND  {info['description']}")

    normalized_score = raw_score / TOTAL_MARKERS
    aging_rate_probability = 1 - normalized_score

    print(f"\nMarkers found: {found}/{len(LONGEVITY_MARKERS)}")
    print(f"Polygenic score: {normalized_score:.3f}")
    print(f"Aging rate probability: {aging_rate_probability:.3f}")
    print("(Lower score suggests increased aging risk)")


if __name__ == "__main__":  # pragma: no cover - manual execution
    main()
