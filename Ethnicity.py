import pandas as pd

# Load your genotype data into a DataFrame named `df`
df = pd.read_csv('Genome.txt', sep='\t', comment='#', header=None)
df.columns = ['rsid', 'chromosome', 'position', 'genotype']
# Ethnicity-informative SNPs – variants known to differ across populations
ethnic_snps = {
    'rs1426654': 'SLC24A5 – Skin pigmentation, differentiates European vs. African ancestry',
    'rs16891982': 'SLC45A2 – Skin, hair, eye color variation',
    'rs3827760': 'EDAR – Hair thickness, East Asian vs. other populations',
    'rs12913832': 'HERC2/OCA2 – Blue vs. brown eye color, European ancestry',
    'rs2736100': 'TERT – Telomere length differences across populations',
    'rs6058017': 'HERC2 – European light hair and eye pigmentation',
    'rs11076174': 'SLC39A8 – Height and metabolic traits differing by ancestry',
    'rs334': 'HBB – Sickle cell trait (African ancestry indicator)',
    # add more ancestry-informative markers here
}

print("\n Ethnicity-informative SNPs:\n")
for rsid, note in ethnic_snps.items():
    row = df[df["rsid"] == rsid]
    if not row.empty:
        genotype = row.iloc[0]["genotype"]
        print(f"{rsid}: {genotype} → {note}")
    else:
        print(f"{rsid}: SNP not found in your file")