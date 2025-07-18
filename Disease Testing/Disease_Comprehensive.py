"""
Comprehensive Disease Risk Analysis
Combines detailed disease risk report, statistics, and polygenic risk scores using Genome.txt
Includes population grouping, risk stratification, and interpretation.
"""

import pandas as pd

def load_data():
    df = pd.read_csv("Genome.txt", sep="\t", comment="#", header=None)
    df.columns = ["rsid", "chromosome", "position", "genotype"]
    return df

def get_allele_count(genotype, risk_allele):
    return genotype.count(risk_allele)

def get_disease_stats():
    # Disease SNPs with extra statistics (from Disease_Statistics.py)
    return {
        'rs429358': {'gene': 'APOE', 'population_frequency': 0.15, 'odds_ratio': 3.5, 'confidence_interval': '2.8-4.4'},
        'rs2187668': {'gene': 'HLA-DQ2', 'population_frequency': 0.08, 'odds_ratio': 2.8, 'confidence_interval': '2.1-3.7'},
        'rs4988235': {'gene': 'LCT', 'population_frequency': 0.65, 'odds_ratio': 4.2, 'confidence_interval': '3.5-5.1'},
        'rs1800562': {'gene': 'HFE', 'population_frequency': 0.06, 'odds_ratio': 8.5, 'confidence_interval': '6.2-11.8'},
        'rs6025': {'gene': 'F5', 'population_frequency': 0.05, 'odds_ratio': 3.8, 'confidence_interval': '2.9-5.0'},
        'rs7903146': {'gene': 'TCF7L2', 'population_frequency': 0.30, 'odds_ratio': 1.4, 'confidence_interval': '1.2-1.6'},
        'rs1333049': {'gene': '9p21', 'population_frequency': 0.45, 'odds_ratio': 1.3, 'confidence_interval': '1.1-1.5'},
        'rs10490924': {'gene': 'ARMS2', 'population_frequency': 0.25, 'odds_ratio': 2.1, 'confidence_interval': '1.7-2.6'},
        'rs2066847': {'gene': 'NOD2', 'population_frequency': 0.10, 'odds_ratio': 2.5, 'confidence_interval': '1.9-3.3'},
        'rs2476601': {'gene': 'PTPN22', 'population_frequency': 0.12, 'odds_ratio': 1.8, 'confidence_interval': '1.4-2.3'},
        'rs3135388': {'gene': 'HLA-DRB1', 'population_frequency': 0.15, 'odds_ratio': 2.2, 'confidence_interval': '1.7-2.8'},
        'rs9939609': {'gene': 'FTO', 'population_frequency': 0.40, 'odds_ratio': 1.2, 'confidence_interval': '1.1-1.3'},
        'rs356219': {'gene': 'SNCA', 'population_frequency': 0.20, 'odds_ratio': 1.4, 'confidence_interval': '1.1-1.7'},
        'rs10484554': {'gene': 'HLA-C', 'population_frequency': 0.18, 'odds_ratio': 1.6, 'confidence_interval': '1.3-2.0'},
        'rs10993994': {'gene': 'MSMB', 'population_frequency': 0.35, 'odds_ratio': 1.3, 'confidence_interval': '1.1-1.5'},
        'rs1799950': {'gene': 'BRCA1', 'population_frequency': 0.02, 'odds_ratio': 2.8, 'confidence_interval': '1.8-4.3'},
        'rs6265': {'gene': 'BDNF', 'population_frequency': 0.25, 'odds_ratio': 1.2, 'confidence_interval': '1.0-1.4'},
        'rs12134493': {'gene': 'CACNA1A', 'population_frequency': 0.30, 'odds_ratio': 1.3, 'confidence_interval': '1.1-1.5'},
        'rs7216389': {'gene': '17q21', 'population_frequency': 0.22, 'odds_ratio': 1.4, 'confidence_interval': '1.2-1.7'},
        'rs7574865': {'gene': 'STAT4', 'population_frequency': 0.15, 'odds_ratio': 1.6, 'confidence_interval': '1.3-2.0'},
    }

def disease_risk_report(df):
    # Comprehensive SNP panel (from Comprehensive_Disease.py)
    comprehensive_snps = {
        # Global/common disease markers
        'rs429358':    {'disease': "Alzheimer's (late onset)", 'risk': ['C'], 'description': "APOE ε4 allele increases risk of late-onset Alzheimer's.", 'population': 'Global'},
        'rs2187668':   {'disease': "Celiac disease", 'risk': ['T'], 'description': "Associated with HLA-DQ2, key in immune response to gluten.", 'population': 'Global'},
        'rs4988235':   {'disease': "Lactose intolerance", 'risk': ['C'], 'description': "CC genotype likely causes lactose intolerance.", 'population': 'Global'},
        'rs1800562':   {'disease': "Hemochromatosis", 'risk': ['G'], 'description': "Mutation in HFE gene leads to iron overload.", 'population': 'Global'},
        'rs6025':      {'disease': "Factor V Leiden (thrombophilia)", 'risk': ['A'], 'description': "Increased risk of blood clots.", 'population': 'Global'},
        'rs7903146':   {'disease': "Type 2 Diabetes", 'risk': ['T'], 'description': "TCF7L2 gene variant raises diabetes risk.", 'population': 'Global'},
        'rs1333049':   {'disease': "Coronary artery disease", 'risk': ['C'], 'description': "Strong association with heart disease.", 'population': 'Global'},
        'rs10490924':  {'disease': "Macular degeneration", 'risk': ['T'], 'description': "Increases risk of age-related vision loss.", 'population': 'Global'},
        'rs2066847':   {'disease': "Crohn's disease", 'risk': ['C'], 'description': "Mutation in NOD2 gene linked to Crohn's.", 'population': 'Global'},
        'rs2476601':   {'disease': "Rheumatoid arthritis", 'risk': ['A'], 'description': "PTPN22 gene variant increases autoimmune risk.", 'population': 'Global'},
        'rs3135388':   {'disease': "Multiple sclerosis", 'risk': ['T'], 'description': "HLA-DRB1*15:01 allele linked to MS.", 'population': 'Global'},
        'rs9939609':   {'disease': "Obesity", 'risk': ['A'], 'description': "FTO gene variant associated with higher BMI.", 'population': 'Global'},
        'rs356219':    {'disease': "Parkinson's disease", 'risk': ['G'], 'description': "SNCA gene variant increases risk.", 'population': 'Global'},
        'rs10484554':  {'disease': "Psoriasis", 'risk': ['T'], 'description': "Linked to immune skin response.", 'population': 'Global'},
        'rs10993994':  {'disease': "Prostate cancer", 'risk': ['T'], 'description': "Risk allele in MSMB gene region.", 'population': 'Global'},
        'rs1799950':   {'disease': "Breast cancer (BRCA1 proxy)", 'risk': ['G'], 'description': "Rare variant possibly linked to BRCA1.", 'population': 'Global'},
        'rs6265':      {'disease': "Depression / neuroticism", 'risk': ['C'], 'description': "BDNF gene variant may affect mood regulation.", 'population': 'Global'},
        'rs12134493':  {'disease': "Migraine", 'risk': ['A'], 'description': "Variant in CACNA1A gene affects migraine susceptibility.", 'population': 'Global'},
        'rs7216389':   {'disease': "Asthma (childhood)", 'risk': ['T'], 'description': "Variant on 17q21 linked to early asthma.", 'population': 'Global'},
        'rs7574865':   {'disease': "Lupus (SLE)", 'risk': ['T'], 'description': "STAT4 gene variant common in autoimmunity.", 'population': 'Global'},
        # East Asian-specific and high-prevalence markers
        'rs671':       {'disease': "Alcohol flush reaction (ALDH2 deficiency)", 'risk': ['A'], 'description': "ALDH2*2 allele causes alcohol intolerance, common in East Asians.", 'population': 'East Asian'},
        'rs1229984':   {'disease': "Alcohol metabolism (ADH1B)", 'risk': ['A'], 'description': "ADH1B*2 allele increases alcohol metabolism, common in East Asians.", 'population': 'East Asian'},
        'rs2075650':   {'disease': "Alzheimer's disease (APOE region, East Asian)", 'risk': ['G'], 'description': "Associated with Alzheimer's in East Asians.", 'population': 'East Asian'},
        'rs1801133':   {'disease': "Homocysteine metabolism (MTHFR)", 'risk': ['T'], 'description': "MTHFR C677T variant, higher risk of hyperhomocysteinemia, common in East Asians.", 'population': 'East Asian'},
        'rs2231142':   {'disease': "Gout (ABCG2)", 'risk': ['T'], 'description': "ABCG2 Q141K variant, high gout risk in East Asians.", 'population': 'East Asian'},
        'rs2285666':   {'disease': "ACE2 expression (COVID-19 susceptibility)", 'risk': ['A'], 'description': "Variant may affect ACE2 expression, studied in East Asians.", 'population': 'East Asian'},
        'rs11200638':  {'disease': "Age-related macular degeneration (HTRA1)", 'risk': ['A'], 'description': "HTRA1 risk allele, high prevalence in East Asians.", 'population': 'East Asian'},
        'rs1800414':   {'disease': "Skin pigmentation (OCA2)", 'risk': ['G'], 'description': "OCA2 variant, common in East Asians.", 'population': 'East Asian'},
        'rs1042522':   {'disease': "Cancer risk (TP53)", 'risk': ['C'], 'description': "TP53 Arg72Pro, cancer risk variant, higher in East Asians.", 'population': 'East Asian'},
        'rs11655237':  {'disease': "Gastric cancer (LINC00673)", 'risk': ['A'], 'description': "LINC00673 variant, gastric cancer risk in East Asians.", 'population': 'East Asian'},
        'rs2736100':   {'disease': "Lung cancer (TERT)", 'risk': ['A'], 'description': "TERT variant, lung cancer risk in East Asians.", 'population': 'East Asian'},
        'rs1801274':   {'disease': "Autoimmune disease (FCGR2A)", 'risk': ['A'], 'description': "FCGR2A variant, SLE/autoimmunity risk, higher in East Asians.", 'population': 'East Asian'},
    }
    stats = get_disease_stats()
    print("\n==============================")
    print("COMPREHENSIVE DISEASE RISK ANALYSIS")
    print("==============================\n")
    populations = {}
    for rsid, data in comprehensive_snps.items():
        pop = data['population']
        if pop not in populations:
            populations[pop] = []
        populations[pop].append((rsid, data))
    total_snps = len(comprehensive_snps)
    found_snps = 0
    risk_alleles_found = 0
    for population in ['Global', 'East Asian']:
        if population in populations:
            print(f"{'='*60}")
            print(f"POPULATION: {population}")
            print(f"{'='*60}\n")
            for rsid, data in populations[population]:
                match = df[df['rsid'] == rsid]
                if not match.empty:
                    genotype = match.iloc[0]['genotype']
                    has_risk = any(allele in genotype for allele in data['risk'])
                    risk_status = "RISK DETECTED" if has_risk else "NO RISK"
                    found_snps += 1
                    if has_risk:
                        risk_alleles_found += 1
                    print(f"{'-'*50}")
                    print(f"Disease: {data['disease']}")
                    print(f"{'-'*50}")
                    print(f"SNP: {rsid}")
                    print(f"Your Genotype: {genotype}")
                    print(f"Risk Allele(s): {', '.join(data['risk'])}")
                    print(f"Status: {risk_status}")
                    print(f"Description: {data['description']}")
                    # Add extra stats if available
                    if rsid in stats:
                        s = stats[rsid]
                        print(f"Gene: {s['gene']}")
                        print(f"Population Frequency: {s['population_frequency']}")
                        print(f"Odds Ratio: {s['odds_ratio']}")
                        print(f"Confidence Interval: {s['confidence_interval']}")
                    print()
                else:
                    print(f"{'-'*50}")
                    print(f"Disease: {data['disease']}")
                    print(f"{'-'*50}")
                    print(f"SNP: {rsid} - NOT FOUND IN DATA")
                    print(f"Description: {data['description']}")
                    if rsid in stats:
                        s = stats[rsid]
                        print(f"Gene: {s['gene']}")
                        print(f"Population Frequency: {s['population_frequency']}")
                        print(f"Odds Ratio: {s['odds_ratio']}")
                        print(f"Confidence Interval: {s['confidence_interval']}")
                    print()
    print(f"\nSUMMARY STATISTICS")
    print(f"{'='*60}")
    print(f"Total SNPs analyzed: {total_snps}")
    print(f"SNPs found in your data: {found_snps}")
    print(f"Risk alleles found: {risk_alleles_found}")
    print(f"Data coverage: {(found_snps/total_snps)*100:.1f}%")

def prs_analysis(df):
    # PRS models (from Comprehensive_Disease.py)
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
        'Alzheimer\'s Disease': {
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

def interpretation():
    print("==============================")
    print("INTERPRETATION")
    print("==============================")
    print("- RISK DETECTED: You carry risk alleles for this condition")
    print("- NO RISK: You do not carry the known risk alleles")
    print("- NOT FOUND: This SNP was not present in your genetic data")
    print("- PRS: Higher scores = higher genetic risk compared to population average")
    print("- These results are for educational purposes only")
    print("- Consult healthcare professionals for medical advice")

def main():
    df = load_data()
    disease_risk_report(df)
    prs_analysis(df)
    interpretation()

if __name__ == "__main__":
    main()