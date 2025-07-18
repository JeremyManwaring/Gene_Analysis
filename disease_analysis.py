import pandas as pd
import random

# Load genotype data

def load_genome_dataframe(filename="Genome.txt"):
    """Load genome data into a pandas DataFrame."""
    try:
        df = pd.read_csv(filename, sep='\t', comment='#', header=None)
        df.columns = ['rsid', 'chromosome', 'position', 'genotype']
    except FileNotFoundError:
        print(f"Warning: {filename} not found. Returning empty DataFrame.")
        df = pd.DataFrame(columns=['rsid', 'chromosome', 'position', 'genotype'])
    return df


def load_genome_dict(filename="Genome.txt", snp_list=None):
    """Load genome data into a dictionary or create simulated data."""
    genome_data = {}
    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 2:
                        snp_id = parts[0]
                        genotype = parts[1]
                        genome_data[snp_id] = genotype
    except FileNotFoundError:
        print(f"Warning: {filename} not found. Using simulated data.")
        if snp_list is None:
            snp_list = []
        genotypes = ['AA', 'AG', 'GG', 'AC', 'CC', 'AT', 'TT', 'CT', 'TC']
        for snp in snp_list:
            genome_data[snp] = random.choice(genotypes)
    return genome_data


def get_allele_count(genotype, risk_allele):
    """Count the number of risk alleles in a genotype string."""
    return genotype.count(risk_allele)

# Disease SNPs (combined from previous modules)
comprehensive_snps = {
    # Global/common disease markers
    'rs429358':    {'disease': "Alzheimer's (late onset)", 'risk': ['C'], 'description': "APOE ε4 allele increases risk of late-onset Alzheimer's."},
    'rs2187668':   {'disease': "Celiac disease", 'risk': ['T'], 'description': "Associated with HLA-DQ2, key in immune response to gluten."},
    'rs4988235':   {'disease': "Lactose intolerance", 'risk': ['C'], 'description': "CC genotype likely causes lactose intolerance."},
    'rs1800562':   {'disease': "Hemochromatosis", 'risk': ['G'], 'description': "Mutation in HFE gene leads to iron overload."},
    'rs6025':      {'disease': "Factor V Leiden (thrombophilia)", 'risk': ['A'], 'description': "Increased risk of blood clots."},
    'rs7903146':   {'disease': "Type 2 Diabetes", 'risk': ['T'], 'description': "TCF7L2 gene variant raises diabetes risk."},
    'rs1333049':   {'disease': "Coronary artery disease", 'risk': ['C'], 'description': "Strong association with heart disease."},
    'rs10490924':  {'disease': "Macular degeneration", 'risk': ['T'], 'description': "Increases risk of age-related vision loss."},
    'rs2066847':   {'disease': "Crohn's disease", 'risk': ['C'], 'description': "Mutation in NOD2 gene linked to Crohn's."},
    'rs2476601':   {'disease': "Rheumatoid arthritis", 'risk': ['A'], 'description': "PTPN22 gene variant increases autoimmune risk."},
    'rs3135388':   {'disease': "Multiple sclerosis", 'risk': ['T'], 'description': "HLA-DRB1*15:01 allele linked to MS."},
    'rs9939609':   {'disease': "Obesity", 'risk': ['A'], 'description': "FTO gene variant associated with higher BMI."},
    'rs356219':    {'disease': "Parkinson's disease", 'risk': ['G'], 'description': "SNCA gene variant increases risk."},
    'rs10484554':  {'disease': "Psoriasis", 'risk': ['T'], 'description': "Linked to immune skin response."},
    'rs10993994':  {'disease': "Prostate cancer", 'risk': ['T'], 'description': "Risk allele in MSMB gene region."},
    'rs1799950':   {'disease': "Breast cancer (BRCA1 proxy)", 'risk': ['G'], 'description': "Rare variant possibly linked to BRCA1."},
    'rs6265':      {'disease': "Depression / neuroticism", 'risk': ['C'], 'description': "BDNF gene variant may affect mood regulation."},
    'rs12134493':  {'disease': "Migraine", 'risk': ['A'], 'description': "Variant in CACNA1A gene affects migraine susceptibility."},
    'rs7216389':   {'disease': "Asthma (childhood)", 'risk': ['T'], 'description': "Variant on 17q21 linked to early asthma."},
    'rs7574865':   {'disease': "Lupus (SLE)", 'risk': ['T'], 'description': "STAT4 gene variant common in autoimmunity."},
    # East Asian-specific markers
    'rs671':       {'disease': "Alcohol flush reaction (ALDH2 deficiency)", 'risk': ['A'], 'description': "ALDH2*2 allele causes alcohol intolerance, common in East Asians."},
    'rs1229984':   {'disease': "Alcohol metabolism (ADH1B)", 'risk': ['A'], 'description': "ADH1B*2 allele increases alcohol metabolism, common in East Asians."},
    'rs2075650':   {'disease': "Alzheimer's disease (APOE region, East Asian)", 'risk': ['G'], 'description': "Associated with Alzheimer's in East Asians."},
    'rs1801133':   {'disease': "Homocysteine metabolism (MTHFR)", 'risk': ['T'], 'description': "MTHFR C677T variant, higher risk of hyperhomocysteinemia, common in East Asians."},
    'rs2231142':   {'disease': "Gout (ABCG2)", 'risk': ['T'], 'description': "ABCG2 Q141K variant, high gout risk in East Asians."},
    'rs2285666':   {'disease': "ACE2 expression (COVID-19 susceptibility)", 'risk': ['A'], 'description': "Variant may affect ACE2 expression, studied in East Asians."},
    'rs11200638':  {'disease': "Age-related macular degeneration (HTRA1)", 'risk': ['A'], 'description': "HTRA1 risk allele, high prevalence in East Asians."},
    'rs1800414':   {'disease': "Skin pigmentation (OCA2)", 'risk': ['G'], 'description': "OCA2 variant, common in East Asians."},
    'rs1042522':   {'disease': "Cancer risk (TP53)", 'risk': ['C'], 'description': "TP53 Arg72Pro, cancer risk variant, higher in East Asians."},
    'rs11655237':  {'disease': "Gastric cancer (LINC00673)", 'risk': ['A'], 'description': "LINC00673 variant, gastric cancer risk in East Asians."},
    'rs2736100':   {'disease': "Lung cancer (TERT)", 'risk': ['A'], 'description': "TERT variant, lung cancer risk in East Asians."},
    'rs1801274':   {'disease': "Autoimmune disease (FCGR2A)", 'risk': ['A'], 'description': "FCGR2A variant, SLE/autoimmunity risk, higher in East Asians."},
}

# PRS models for major diseases
prs_models = {
    'Coronary Artery Disease': {
        'rs1333049': ('C', 0.25, '9p21 locus'),
        'rs10757278': ('G', 0.18, '9p21 locus'),
        'rs2383206': ('G', 0.15, '9p21 locus'),
        'rs2383207': ('A', 0.12, '9p21 locus'),
        'rs10757274': ('G', 0.10, '9p21 locus'),
    },
    'Type 2 Diabetes': {
        'rs7903146': ('T', 0.30, 'TCF7L2'),
        'rs1801282': ('G', 0.20, 'PPARG'),
        'rs5219': ('T', 0.15, 'KCNJ11'),
        'rs13266634': ('C', 0.12, 'SLC30A8'),
        'rs4402960': ('T', 0.10, 'IGF2BP2'),
    },
    "Alzheimer's Disease": {
        'rs429358': ('C', 0.40, 'APOE'),
        'rs7412': ('C', 0.35, 'APOE'),
        'rs2075650': ('G', 0.25, 'TOMM40'),
        'rs157580': ('G', 0.20, 'TOMM40'),
        'rs11556505': ('T', 0.18, 'PICALM'),
    },
    'Breast Cancer': {
        'rs2981582': ('C', 0.25, 'FGFR2'),
        'rs3803662': ('C', 0.20, 'TNRC9'),
        'rs889312': ('C', 0.18, 'MAP3K1'),
        'rs3817198': ('T', 0.16, 'LSP1'),
        'rs13281615': ('G', 0.14, '8q24'),
    },
    'Prostate Cancer': {
        'rs10993994': ('T', 0.30, 'MSMB'),
        'rs7931342': ('T', 0.25, '11q13'),
        'rs2735839': ('G', 0.20, 'KLK3'),
        'rs17632542': ('T', 0.18, 'KLK3'),
        'rs1859962': ('G', 0.16, '17q24.3'),
    },
    'Depression': {
        'rs6265': ('C', 0.25, 'BDNF'),
        'rs1360780': ('T', 0.20, 'FKBP5'),
        'rs3800373': ('C', 0.18, 'FKBP5'),
        'rs9470080': ('T', 0.16, 'SLC6A4'),
        'rs25531': ('A', 0.14, 'SLC6A4'),
    },
}

# Detailed SNP info for statistics
# (using the dictionary from Disease_Statistics)
detailed_snps = {
    'rs429358': {
        'disease': "Alzheimer's (late onset)",
        'gene': 'APOE',
        'risk_alleles': ['C'],
        'risk_description': "APOE ε4 allele increases risk of late-onset Alzheimer's",
        'population_frequency': 0.15,
        'odds_ratio': 3.5,
        'confidence_interval': '2.8-4.4'
    },
    'rs2187668': {
        'disease': "Celiac disease",
        'gene': 'HLA-DQ2',
        'risk_alleles': ['T'],
        'risk_description': "Associated with HLA-DQ2, key in immune response to gluten",
        'population_frequency': 0.08,
        'odds_ratio': 2.8,
        'confidence_interval': '2.1-3.7'
    },
    'rs4988235': {
        'disease': "Lactose intolerance",
        'gene': 'LCT',
        'risk_alleles': ['C'],
        'risk_description': "CC genotype likely causes lactose intolerance",
        'population_frequency': 0.65,
        'odds_ratio': 4.2,
        'confidence_interval': '3.5-5.1'
    },
    'rs1800562': {
        'disease': "Hemochromatosis",
        'gene': 'HFE',
        'risk_alleles': ['G'],
        'risk_description': "Mutation in HFE gene leads to iron overload",
        'population_frequency': 0.06,
        'odds_ratio': 8.5,
        'confidence_interval': '6.2-11.8'
    },
    'rs6025': {
        'disease': "Factor V Leiden (thrombophilia)",
        'gene': 'F5',
        'risk_alleles': ['A'],
        'risk_description': "Increased risk of blood clots",
        'population_frequency': 0.05,
        'odds_ratio': 3.8,
        'confidence_interval': '2.9-5.0'
    },
    'rs7903146': {
        'disease': "Type 2 Diabetes",
        'gene': 'TCF7L2',
        'risk_alleles': ['T'],
        'risk_description': "TCF7L2 gene variant raises diabetes risk",
        'population_frequency': 0.30,
        'odds_ratio': 1.4,
        'confidence_interval': '1.2-1.6'
    },
    'rs1333049': {
        'disease': "Coronary artery disease",
        'gene': '9p21',
        'risk_alleles': ['C'],
        'risk_description': "Strong association with heart disease",
        'population_frequency': 0.45,
        'odds_ratio': 1.3,
        'confidence_interval': '1.1-1.5'
    },
    'rs10490924': {
        'disease': "Macular degeneration",
        'gene': 'ARMS2',
        'risk_alleles': ['T'],
        'risk_description': "Increases risk of age-related vision loss",
        'population_frequency': 0.25,
        'odds_ratio': 2.1,
        'confidence_interval': '1.7-2.6'
    },
    'rs2066847': {
        'disease': "Crohn's disease",
        'gene': 'NOD2',
        'risk_alleles': ['C'],
        'risk_description': "Mutation in NOD2 gene linked to Crohn's",
        'population_frequency': 0.10,
        'odds_ratio': 2.5,
        'confidence_interval': '1.9-3.3'
    },
    'rs2476601': {
        'disease': "Rheumatoid arthritis",
        'gene': 'PTPN22',
        'risk_alleles': ['A'],
        'risk_description': "PTPN22 gene variant increases autoimmune risk",
        'population_frequency': 0.12,
        'odds_ratio': 1.8,
        'confidence_interval': '1.4-2.3'
    },
    'rs3135388': {
        'disease': "Multiple sclerosis",
        'gene': 'HLA-DRB1',
        'risk_alleles': ['T'],
        'risk_description': "HLA-DRB1*15:01 allele linked to MS",
        'population_frequency': 0.15,
        'odds_ratio': 2.2,
        'confidence_interval': '1.7-2.8'
    },
    'rs9939609': {
        'disease': "Obesity",
        'gene': 'FTO',
        'risk_alleles': ['A'],
        'risk_description': "FTO gene variant associated with higher BMI",
        'population_frequency': 0.40,
        'odds_ratio': 1.2,
        'confidence_interval': '1.1-1.3'
    },
    'rs356219': {
        'disease': "Parkinson's disease",
        'gene': 'SNCA',
        'risk_alleles': ['G'],
        'risk_description': "SNCA gene variant increases risk",
        'population_frequency': 0.20,
        'odds_ratio': 1.4,
        'confidence_interval': '1.1-1.7'
    },
    'rs10484554': {
        'disease': "Psoriasis",
        'gene': 'HLA-C',
        'risk_alleles': ['T'],
        'risk_description': "Linked to immune skin response",
        'population_frequency': 0.18,
        'odds_ratio': 1.6,
        'confidence_interval': '1.3-2.0'
    },
    'rs10993994': {
        'disease': "Prostate cancer",
        'gene': 'MSMB',
        'risk_alleles': ['T'],
        'risk_description': "Risk allele in MSMB gene region",
        'population_frequency': 0.35,
        'odds_ratio': 1.3,
        'confidence_interval': '1.1-1.5'
    },
    'rs1799950': {
        'disease': "Breast cancer (BRCA1 proxy)",
        'gene': 'BRCA1',
        'risk_alleles': ['G'],
        'risk_description': "Rare variant possibly linked to BRCA1",
        'population_frequency': 0.02,
        'odds_ratio': 2.8,
        'confidence_interval': '1.8-4.3'
    },
    'rs6265': {
        'disease': "Depression / neuroticism",
        'gene': 'BDNF',
        'risk_alleles': ['C'],
        'risk_description': "BDNF gene variant may affect mood regulation",
        'population_frequency': 0.25,
        'odds_ratio': 1.2,
        'confidence_interval': '1.0-1.4'
    },
    'rs12134493': {
        'disease': "Migraine",
        'gene': 'CACNA1A',
        'risk_alleles': ['A'],
        'risk_description': "Variant in CACNA1A gene affects migraine susceptibility",
        'population_frequency': 0.30,
        'odds_ratio': 1.3,
        'confidence_interval': '1.1-1.5'
    },
    'rs7216389': {
        'disease': "Asthma (childhood)",
        'gene': '17q21',
        'risk_alleles': ['T'],
        'risk_description': "Variant on 17q21 linked to early asthma",
        'population_frequency': 0.22,
        'odds_ratio': 1.4,
        'confidence_interval': '1.2-1.7'
    },
    'rs7574865': {
        'disease': "Lupus (SLE)",
        'gene': 'STAT4',
        'risk_alleles': ['T'],
        'risk_description': "STAT4 gene variant common in autoimmunity",
        'population_frequency': 0.15,
        'odds_ratio': 1.6,
        'confidence_interval': '1.3-2.0'
    },
}


def risk_summary(filename="Genome.txt"):
    """Print a simple risk summary for known disease SNPs."""
    df = load_genome_dataframe(filename)
    results = []
    print("\n GENETIC RISK ANALYSIS:\n")
    for rsid, data in comprehensive_snps.items():
        match = df[df['rsid'] == rsid]
        if not match.empty:
            genotype = match.iloc[0]['genotype']
            has_risk = any(allele in genotype for allele in data['risk'])
            results.append({
                'disease': data['disease'],
                'snp': rsid,
                'genotype': genotype,
                'risk_alleles': data['risk'],
                'has_risk': has_risk,
                'description': data['description']
            })
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
            results.append({
                'disease': data['disease'],
                'snp': rsid,
                'genotype': None,
                'risk_alleles': data['risk'],
                'has_risk': False,
                'description': data['description']
            })
            print(f"Disease: {data['disease']}")
            print(f"  - SNP: {rsid} not found in your raw data.\n")
    return results


def disease_statistics(filename="Genome.txt"):
    """Analyze disease risk and print detailed statistics."""
    genome_data = load_genome_dict(filename, detailed_snps.keys())

    print("📊 DISEASE RISK STATISTICS")
    print("=" * 60)

    total_snps = len(detailed_snps)
    found_snps = 0
    risk_alleles_found = 0
    high_risk_diseases = []
    moderate_risk_diseases = []
    low_risk_diseases = []

    print(f"\n🔬 ANALYZING {total_snps} DISEASE-ASSOCIATED SNPS")
    print("-" * 60)

    for snp_id, data in detailed_snps.items():
        if snp_id in genome_data:
            found_snps += 1
            genotype = genome_data[snp_id]

            has_risk_allele = any(allele in genotype for allele in data['risk_alleles'])
            if has_risk_allele:
                risk_alleles_found += 1

            odds_ratio = data['odds_ratio']
            if odds_ratio >= 3.0:
                risk_level = "HIGH"
                high_risk_diseases.append(data['disease'])
            elif odds_ratio >= 1.5:
                risk_level = "MODERATE"
                moderate_risk_diseases.append(data['disease'])
            else:
                risk_level = "LOW"
                low_risk_diseases.append(data['disease'])

            print(f"\n🏥 {data['disease'].upper()}")
            print(f"   SNP: {snp_id}")
            print(f"   Gene: {data['gene']}")
            print(f"   Your Genotype: {genotype}")
            print(f"   Risk Alleles: {', '.join(data['risk_alleles'])}")
            print(f"   Population Frequency: {data['population_frequency']:.1%}")
            print(f"   Odds Ratio: {data['odds_ratio']} (95% CI: {data['confidence_interval']})")
            print(f"   Risk Level: {risk_level}")
            print(f"   Risk Allele Present: {'YES' if has_risk_allele else 'NO'}")
            print(f"   Description: {data['risk_description']}")
        else:
            print(f"\n🏥 {data['disease'].upper()}")
            print(f"   SNP: {snp_id} - NOT FOUND IN GENOME DATA")

    print(f"\n📈 SUMMARY STATISTICS")
    print("=" * 60)
    print(f"Total SNPs Analyzed: {total_snps}")
    print(f"SNPs Found in Genome: {found_snps} ({found_snps/total_snps:.1%})")
    print(f"Risk Alleles Detected: {risk_alleles_found}")
    print(f"Risk Detection Rate: {risk_alleles_found/found_snps:.1%}" if found_snps > 0 else "Risk Detection Rate: N/A")

    print(f"\n🎯 RISK CATEGORIES")
    print(f"High Risk Diseases: {len(high_risk_diseases)}")
    for disease in high_risk_diseases:
        print(f"  • {disease}")

    print(f"\nModerate Risk Diseases: {len(moderate_risk_diseases)}")
    for disease in moderate_risk_diseases:
        print(f"  • {disease}")

    print(f"\nLow Risk Diseases: {len(low_risk_diseases)}")
    for disease in low_risk_diseases:
        print(f"  • {disease}")

    print(f"\n📊 POPULATION COMPARISON")
    print("=" * 60)
    print("Your genetic profile compared to general population:")
    print(f"• SNPs with risk alleles: {risk_alleles_found}/{found_snps}")
    high_risk_count = len([d for d in high_risk_diseases if d in [snp_data['disease'] for snp_id, snp_data in detailed_snps.items() if snp_id in genome_data and any(allele in genome_data[snp_id] for allele in snp_data['risk_alleles'])]])
    print(f"• High-risk variants: {high_risk_count}")
    overall = 'Elevated' if risk_alleles_found > found_snps/2 else 'Average' if risk_alleles_found > found_snps/4 else 'Lower'
    print(f"• Overall risk profile: {overall}")

    return {
        'total_snps': total_snps,
        'found_snps': found_snps,
        'risk_alleles_found': risk_alleles_found,
        'high_risk_diseases': high_risk_diseases,
        'moderate_risk_diseases': moderate_risk_diseases,
        'low_risk_diseases': low_risk_diseases
    }


def polygenic_risk_scores(filename="Genome.txt"):
    """Calculate polygenic risk scores for supported traits."""
    df = load_genome_dataframe(filename)
    results = {}

    print(f"{'='*60}")
    print("POLYGENIC RISK SCORES (PRS)")
    print(f"{'='*60}\n")

    for trait, snps in prs_models.items():
        print(f"{'='*50}")
        print(f"TRAIT: {trait}")
        print(f"{'='*50}")

        prs = 0.0
        snps_found = 0
        snp_details = []

        for rsid, (risk_allele, weight, source) in snps.items():
            match = df[df['rsid'] == rsid]
            if not match.empty:
                genotype = match.iloc[0]['genotype']
                count = get_allele_count(genotype, risk_allele)
                prs += count * weight
                snps_found += 1
                snp_details.append((rsid, genotype, risk_allele, weight, count, source))
            else:
                snp_details.append((rsid, 'Not found', risk_allele, weight, 0, source))

        print(f"{'SNP':<12}{'Genotype':<12}{'Risk':<8}{'Weight':<10}{'Count':<8}{'Source':<15}")
        print(f"{'-'*70}")
        for rsid, genotype, risk_allele, weight, count, source in snp_details:
            print(f"{rsid:<12}{genotype:<12}{risk_allele:<8}{weight:<10}{count:<8}{source:<15}")
        print(f"{'-'*70}")
        print(f"PRS for {trait}: {prs:.2f} (based on {snps_found}/{len(snps)} SNPs found)")
        print()

        results[trait] = prs

    return results


if __name__ == "__main__":
    risk_summary()
    disease_statistics()
    polygenic_risk_scores()
