"""Quick check of key longevity-related SNPs from ``Genome.txt``.

This script highlights a few major markers linked to aging and lifespan and
calculates a simple longevity probability based on those markers.
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
    """Load ``Genome.txt`` as a DataFrame."""
    df = pd.read_csv("Genome.txt", sep="\t", comment="#", header=None)
    df.columns = ["rsid", "chromosome", "position", "genotype"]
    return df


# Major longevity markers with simple weights
LONGEVITY_SNPS = {
    'rs2736100': {'risk': ['A'], 'weight': 0.15, 'description': 'TERT – Telomere maintenance'},
    'rs7726159': {'risk': ['A'], 'weight': 0.12, 'description': 'TERT – Telomere length'},
    'rs1317082': {'risk': ['A'], 'weight': 0.10, 'description': 'TERC – Telomere regulation'},
    'rs1801133': {'risk': ['A'], 'weight': 0.19, 'description': 'MTHFR – Methylation cycle'},
    'rs1042522': {'risk': ['C'], 'weight': 0.10, 'description': 'TP53 – DNA repair'},
}


def genotype_effect(genotype: str, risk_allele: str) -> float:
    """Calculate the effect of a genotype for a given risk allele."""
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
    print("LONGEVITY MARKER SUMMARY")
    print("==============================\n")

    total_weight = sum(snp['weight'] for snp in LONGEVITY_SNPS.values())
    score = 0.0

    for rsid, info in LONGEVITY_SNPS.items():
        row = df[df['rsid'] == rsid]
        if not row.empty:
            genotype = row.iloc[0]['genotype']
            effect = genotype_effect(genotype, info['risk'][0])
            contribution = effect * info['weight']
            score += contribution
            print(f"SNP: {rsid}  Genotype: {genotype}  Effect: {contribution:.3f}  {info['description']}")
        else:
            print(f"SNP: {rsid} - NOT FOUND  {info['description']}")

    probability = (score / total_weight) * 100 if total_weight else 0
    print(f"\nLongevity probability: {probability:.1f}%")
    print("(Higher percentage suggests genetics associated with longer lifespan)")


if __name__ == "__main__":  # pragma: no cover - manual execution
    main()
