"""Comprehensive ancestry marker check using ``Genome.txt``."""

import pandas as pd

df = pd.read_csv("Genome.txt", sep="\t", comment="#", header=None)
df.columns = ["rsid", "chromosome", "position", "genotype"]

# Comprehensive ancestry-informative SNPs
comprehensive_ethnic_snps = {
    # East Asian specific markers
    'rs3827760': 'EDAR - East Asian hair thickness',
    'rs1426654': 'SLC24A5 - Skin pigmentation',
    'rs16891982': 'SLC45A2 - Skin pigmentation',
    'rs12913832': 'HERC2/OCA2 - Eye color',
    'rs1805007': 'MC1R - Red hair',
    'rs1805008': 'MC1R - Red hair variant',
    'rs885479': 'IRF4 - Hair color',
    'rs12203592': 'IRF4 - Hair/skin color',
    'rs1545397': 'OCA2 - Eye color',
    'rs1800401': 'OCA2 - Eye color variant',
    'rs4778138': 'OCA2 - Eye color',
    'rs7495174': 'OCA2 - Eye color',
    'rs4778241': 'OCA2 - Eye color',

    # South Asian markers
    'rs1042602': 'TYR - Skin pigmentation',
    'rs1800414': 'OCA2 - Eye color',
    'rs1800407': 'OCA2 - Eye color',
    'rs1800408': 'OCA2 - Eye color',
    'rs1800409': 'OCA2 - Eye color',
    'rs1800410': 'OCA2 - Eye color',

    # Additional ancestry markers
    'rs2736100': 'TERT - Telomere length',
    'rs6058017': 'HERC2 - European light pigmentation',
    'rs11076174': 'SLC39A8 - Height/metabolic traits',
    'rs1801133': 'MTHFR - Folate metabolism',
    'rs1801131': 'MTHFR - Folate metabolism',
    'rs1805087': 'MTHFR - Folate metabolism',

    # Lactose tolerance markers
    'rs4988235': 'LCT - Lactose tolerance',
    'rs182549': 'LCT - Lactose tolerance',

    # Alcohol metabolism markers
    'rs1229984': 'ADH1B - Alcohol metabolism',
    'rs671': 'ALDH2 - Alcohol metabolism'
}

print("\n==============================")
print("COMPREHENSIVE ANCESTRY ANALYSIS")
print("==============================\n")

# Count matches by category
east_asian_markers = ['rs3827760']
european_markers = ['rs1426654', 'rs16891982', 'rs12913832', 'rs1805007', 'rs1805008', 'rs885479', 'rs12203592', 'rs6058017']
south_asian_markers = ['rs1042602', 'rs1800414', 'rs1800407', 'rs1800408', 'rs1800409', 'rs1800410']
african_markers = ['rs1426654', 'rs16891982', 'rs12913832']
native_american_markers = ['rs1805007', 'rs1805008', 'rs885479', 'rs12203592']
middle_eastern_markers = ['rs1426654', 'rs16891982', 'rs12913832']

found_markers = []
not_found = []

print(f"{'='*60}")
print("ANCESTRY-INFORMATIVE MARKERS")
print(f"{'='*60}\n")

for rsid, description in comprehensive_ethnic_snps.items():
    row = df[df["rsid"] == rsid]
    if not row.empty:
        genotype = row.iloc[0]["genotype"]
        found_markers.append((rsid, genotype, description))
        print(f"{'-'*50}")
        print(f"SNP: {rsid}")
        print(f"Genotype: {genotype}")
        print(f"Function: {description}")
        print()
    else:
        not_found.append(rsid)

print(f"{'='*60}")
print("SUMMARY STATISTICS")
print(f"{'='*60}")
print(f"Markers found: {len(found_markers)} out of {len(comprehensive_ethnic_snps)}")
print(f"Markers missing: {len(not_found)}")

# Analyze patterns
print(f"\n{'='*60}")
print("POPULATION-SPECIFIC COUNTS")
print(f"{'='*60}")

east_asian_count = sum(1 for rsid, _, _ in found_markers if rsid in east_asian_markers)
european_count = sum(1 for rsid, _, _ in found_markers if rsid in european_markers)
south_asian_count = sum(1 for rsid, _, _ in found_markers if rsid in south_asian_markers)
african_count = sum(1 for rsid, _, _ in found_markers if rsid in african_markers)

print(f"East Asian markers found: {east_asian_count}")
print(f"European markers found: {european_count}")
print(f"South Asian markers found: {south_asian_count}")
print(f"African markers found: {african_count}")

print(f"\n{'='*60}")
print("INTERPRETATION")
print(f"{'='*60}")
print("- This analysis examines genetic markers that differ between populations")
print("- Higher counts in a category suggest genetic similarity to that population")
print("- Remember: Genetic ancestry ≠ Cultural identity")
print("- These results are statistical associations, not definitive ancestry")
print("- Individual genetic variation is complex and doesn't fit simple categories")


if __name__ == "__main__":
    pass