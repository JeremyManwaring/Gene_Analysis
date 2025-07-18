import pandas as pd
df = pd.read_csv('Genome.txt', sep='\t', comment='#', header=None)
df.columns = ['rsid', 'chromosome', 'position', 'genotype']
snps = {
    'rs429358':    {'disease': "Alzheimer's (late onset)", 'risk': ['C'], 'description': "APOE ε4 allele increases risk of late-onset Alzheimer's."},
    'rs2187668':   {'disease': "Celiac disease", 'risk': ['T'], 'description': "Associated with HLA-DQ2, key in immune response to gluten."},
    'rs4988235':   {'disease': "Lactose intolerance", 'risk': ['C'], 'description': "CC genotype likely causes lactose intolerance."},
    'rs1800562':   {'disease': "Hemochromatosis", 'risk': ['G'], 'description': "Mutation in HFE gene leads to iron overload."},
    'rs6025':      {'disease': "Factor V Leiden (thrombophilia)", 'risk': ['A'], 'description': "Increased risk of blood clots."},
    'rs7903146':   {'disease': "Type 2 Diabetes", 'risk': ['T'], 'description': "TCF7L2 gene variant raises diabetes risk."},
    'rs1333049':   {'disease': "Coronary artery disease", 'risk': ['C'], 'description': "Strong association with heart disease."},
    'rs10490924':  {'disease': "Macular degeneration", 'risk': ['T'], 'description': "Increases risk of age-related vision loss."},
    'rs2066847':   {'disease': "Crohn's disease", 'risk': ['C'], 'description': "Mutation in NOD2 gene linked to Crohn’s."},
    'rs2476601':   {'disease': "Rheumatoid arthritis", 'risk': ['A'], 'description': "PTPN22 gene variant increases autoimmune risk."},
    'rs3135388':   {'disease': "Multiple sclerosis", 'risk': ['T'], 'description': "HLA-DRB1*15:01 allele linked to MS."},
    'rs9939609':   {'disease': "Obesity", 'risk': ['A'], 'description': "FTO gene variant associated with higher BMI."},
    'rs356219':    {'disease': "Parkinson’s disease", 'risk': ['G'], 'description': "SNCA gene variant increases risk."},
    'rs10484554':  {'disease': "Psoriasis", 'risk': ['T'], 'description': "Linked to immune skin response."},
    'rs10993994':  {'disease': "Prostate cancer", 'risk': ['T'], 'description': "Risk allele in MSMB gene region."},
    'rs1799950':   {'disease': "Breast cancer (BRCA1 proxy)", 'risk': ['G'], 'description': "Rare variant possibly linked to BRCA1."},
    'rs6265':      {'disease': "Depression / neuroticism", 'risk': ['C'], 'description': "BDNF gene variant may affect mood regulation."},
    'rs12134493':  {'disease': "Migraine", 'risk': ['A'], 'description': "Variant in CACNA1A gene affects migraine susceptibility."},
    'rs7216389':   {'disease': "Asthma (childhood)", 'risk': ['T'], 'description': "Variant on 17q21 linked to early asthma."},
    'rs7574865':   {'disease': "Lupus (SLE)", 'risk': ['T'], 'description': "STAT4 gene variant common in autoimmunity."},
    'rs662799':   {'disease': "Hypertriglyceridemia", 'risk': ['C'], 'description': "APOA5 promoter variant linked to high triglycerides."}
}

# Polygenic risk score SNPs – risk variants for complex traits
prs_snps = {
    'rs429358': 'APOE ε4 allele – Alzheimer’s disease risk',
    'rs9939609': 'FTO – Obesity predisposition',
    'rs7903146': 'TCF7L2 – Type 2 diabetes risk',
    'rs1333049': '9p21 locus – Coronary artery disease risk',
    'rs11591147': 'PCSK9 – LDL cholesterol and CAD risk',
    'rs2981582': 'FGFR2 – Breast cancer risk',
    'rs7597593': 'ZNF804A – Schizophrenia risk',
    'rs1067828': 'SLC6A4 – Depression and SSRI response',
    'rs662799': 'APOA5 – Hypertriglyceridemia risk',
    # add more PRS SNPs relevant to your traits here

}
print("\n GENETIC RISK ANALYSIS:\n")
for rsid, data in snps.items():
    match = df[df['rsid'] == rsid]
    if not match.empty:
        genotype = match.iloc[0]['genotype']
        has_risk = any(allele in genotype for allele in data['risk'])
        print(f"Disease: {data['disease']}")
        print(f"  - SNP: {rsid}")
        print(f"  - Your genotype: {genotype}")
        print(f"  - Risk allele(s): {', '.join(data['risk'])}")
        print(f"  - Description: {data['description']}")
        if has_risk:
            print("  Potential genetic risk detected (risk allele present)\n")
        else:
            print("  No known risk alleles found\n")
    else:
        print(f"Disease: {data['disease']}")
        print(f"  - SNP: {rsid} not found in your raw data.\n")