import pandas as pd

# Load your genotype data
df = pd.read_csv('Genome.txt', sep='\t', comment='#', header=None)
df.columns = ['rsid', 'chromosome', 'position', 'genotype']

# Comprehensive ancestry-informative SNPs organized by category
ancestry_categories = {
    'Skin Pigmentation': {
        'rs1426654': 'SLC24A5 - European vs African skin pigmentation',
        'rs16891982': 'SLC45A2 - Skin, hair, eye color variation',
        'rs1042602': 'TYR - Skin pigmentation variation',
    },
    'Eye Color': {
        'rs12913832': 'HERC2/OCA2 - Blue vs brown eye color (European)',
        'rs1545397': 'OCA2 - Eye color variation',
        'rs1800401': 'OCA2 - Eye color variant',
        'rs4778138': 'OCA2 - Eye color',
        'rs7495174': 'OCA2 - Eye color',
        'rs4778241': 'OCA2 - Eye color',
        'rs1800407': 'OCA2 - Eye color',
    },
    'Hair Color & Texture': {
        'rs1805007': 'MC1R - Red hair (European)',
        'rs1805008': 'MC1R - Red hair variant',
        'rs885479': 'IRF4 - Hair color',
        'rs12203592': 'IRF4 - Hair/skin color',
        'rs3827760': 'EDAR - East Asian hair thickness',
    },
    'Metabolic Traits': {
        'rs2736100': 'TERT - Telomere length differences',
        'rs11076174': 'SLC39A8 - Height and metabolic traits',
        'rs1801133': 'MTHFR - Folate metabolism',
        'rs1801131': 'MTHFR - Folate metabolism',
        'rs1805087': 'MTHFR - Folate metabolism',
    },
    'Lactose Tolerance': {
        'rs4988235': 'LCT - Lactose tolerance',
        'rs182549': 'LCT - Lactose tolerance',
    },
    'Alcohol Metabolism': {
        'rs1229984': 'ADH1B - Alcohol metabolism',
        'rs671': 'ALDH2 - Alcohol metabolism',
    },
    'African Ancestry Markers': {
        'rs334': 'HBB - Sickle cell trait indicator',
    },
    'European Ancestry Markers': {
        'rs6058017': 'HERC2 - European light pigmentation',
    }
}

print("\n==============================")
print("COMPREHENSIVE ANCESTRY SUMMARY")
print("==============================\n")

# Track results for summary
all_results = {}
category_counts = {}

for category, snps in ancestry_categories.items():
    print(f"{'='*60}")
    print(f"CATEGORY: {category}")
    print(f"{'='*60}")
    
    category_results = []
    found_count = 0
    
    for rsid, description in snps.items():
        row = df[df["rsid"] == rsid]
        if not row.empty:
            genotype = row.iloc[0]["genotype"]
            found_count += 1
            category_results.append((rsid, genotype, description))
            all_results[rsid] = (genotype, description, category)
            print(f"{'-'*50}")
            print(f"SNP: {rsid}")
            print(f"Genotype: {genotype}")
            print(f"Function: {description}")
            print()
        else:
            print(f"{'-'*50}")
            print(f"SNP: {rsid} - NOT FOUND IN DATA")
            print(f"Function: {description}")
            print()
    
    category_counts[category] = (found_count, len(snps))

# Summary Statistics
print(f"{'='*60}")
print("SUMMARY STATISTICS")
print(f"{'='*60}")

total_found = sum(count for count, total in category_counts.values())
total_snps = sum(total for count, total in category_counts.values())

print(f"Total SNPs analyzed: {total_snps}")
print(f"SNPs found in your data: {total_found}")
print(f"SNPs missing: {total_snps - total_found}")
print(f"Data coverage: {(total_found/total_snps)*100:.1f}%")

print(f"\n{'='*60}")
print("CATEGORY BREAKDOWN")
print(f"{'='*60}")

for category, (found, total) in category_counts.items():
    coverage = (found/total)*100 if total > 0 else 0
    print(f"{category:<25} {found}/{total} SNPs ({coverage:.1f}% coverage)")

# Population-specific analysis
print(f"\n{'='*60}")
print("POPULATION-SPECIFIC ANALYSIS")
print(f"{'='*60}")

# Define population markers
east_asian_markers = ['rs3827760']
european_markers = ['rs1426654', 'rs16891982', 'rs12913832', 'rs1805007', 'rs1805008', 'rs885479', 'rs12203592', 'rs6058017']
african_markers = ['rs1426654', 'rs16891982', 'rs12913832']

east_asian_count = sum(1 for rsid in all_results.keys() if rsid in east_asian_markers)
european_count = sum(1 for rsid in all_results.keys() if rsid in european_markers)
african_count = sum(1 for rsid in all_results.keys() if rsid in african_markers)

print(f"East Asian markers found: {east_asian_count}")
print(f"European markers found: {european_count}")
print(f"African markers found: {african_count}")

# Key findings summary
print(f"\n{'='*60}")
print("KEY FINDINGS")
print(f"{'='*60}")

# Analyze specific patterns
if 'rs3827760' in all_results:
    genotype = all_results['rs3827760'][0]
    if genotype == 'AA' or genotype == 'AG':
        print("✓ East Asian hair thickness variant detected")
    else:
        print("✗ East Asian hair thickness variant not detected")

if 'rs671' in all_results:
    genotype = all_results['rs671'][0]
    if genotype == 'AA' or genotype == 'AG':
        print("✓ Asian flush (ALDH2 deficiency) variant detected")
    else:
        print("✗ Asian flush variant not detected")

if 'rs4988235' in all_results:
    genotype = all_results['rs4988235'][0]
    if 'C' in genotype:
        print("✓ Lactose intolerance variant detected")
    else:
        print("✗ Lactose intolerance variant not detected")

print(f"\n{'='*60}")
print("INTERPRETATION")
print(f"{'='*60}")
print("- This analysis examines genetic markers that differ between populations")
print("- Higher counts in a category suggest genetic similarity to that population")
print("- Remember: Genetic ancestry ≠ Cultural identity")
print("- These results are statistical associations, not definitive ancestry")
print("- Individual genetic variation is complex and doesn't fit simple categories")
print("- Your cultural identity is valid regardless of genetic markers") 