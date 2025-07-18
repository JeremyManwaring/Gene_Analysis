#!/usr/bin/env python3
"""
Disease Risk Statistics
Shows raw SNP data and statistics for disease risk analysis
"""

import random

def load_genome_data(filename):
    """Load genome data from file"""
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
        # Generate simulated data for demonstration
        all_snps = [
            'rs429358', 'rs2187668', 'rs4988235', 'rs1800562', 'rs6025',
            'rs7903146', 'rs1333049', 'rs10490924', 'rs2066847', 'rs2476601',
            'rs3135388', 'rs9939609', 'rs356219', 'rs10484554', 'rs10993994',
            'rs1799950', 'rs6265', 'rs12134493', 'rs7216389', 'rs7574865',
            'rs11591147', 'rs2981582', 'rs7597593', 'rs1067828'
        ]
        for snp in all_snps:
            genotypes = ['AA', 'AG', 'GG', 'AC', 'CC', 'AT', 'TT', 'CT', 'TC']
            genome_data[snp] = random.choice(genotypes)
    
    return genome_data

def analyze_disease_risk():
    """Analyze disease risk with detailed statistics"""
    
    # Disease risk SNPs with detailed information
    disease_snps = {
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
        'rs662799': {
            'disease': 'Hypertriglyceridemia',
            'gene': 'APOA5',
            'risk_alleles': ['C'],
            'risk_description': 'APOA5 promoter variant linked to high triglycerides',
            'population_frequency': 0.20,
            'odds_ratio': 1.9,
            'confidence_interval': '1.4-2.5'
        }
    }
    
    # Load genome data
    genome_data = load_genome_data("Genome.txt")
    
    print("📊 DISEASE RISK STATISTICS")
    print("=" * 60)
    
    # Statistics counters
    total_snps = len(disease_snps)
    found_snps = 0
    risk_alleles_found = 0
    high_risk_diseases = []
    moderate_risk_diseases = []
    low_risk_diseases = []
    
    print(f"\n🔬 ANALYZING {total_snps} DISEASE-ASSOCIATED SNPS")
    print("-" * 60)
    
    for snp_id, data in disease_snps.items():
        if snp_id in genome_data:
            found_snps += 1
            genotype = genome_data[snp_id]
            
            # Check for risk alleles
            has_risk_allele = any(allele in genotype for allele in data['risk_alleles'])
            if has_risk_allele:
                risk_alleles_found += 1
            
            # Calculate risk level based on odds ratio
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
    
    # Summary statistics
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
    
    # Population comparison
    print(f"\n📊 POPULATION COMPARISON")
    print("=" * 60)
    print("Your genetic profile compared to general population:")
    print(f"• SNPs with risk alleles: {risk_alleles_found}/{found_snps}")
    print(f"• High-risk variants: {len([d for d in high_risk_diseases if d in [snp_data['disease'] for snp_id, snp_data in disease_snps.items() if snp_id in genome_data and any(allele in genome_data[snp_id] for allele in snp_data['risk_alleles'])]])}")
    print(f"• Overall risk profile: {'Elevated' if risk_alleles_found > found_snps/2 else 'Average' if risk_alleles_found > found_snps/4 else 'Lower'}")
    
    return {
        'total_snps': total_snps,
        'found_snps': found_snps,
        'risk_alleles_found': risk_alleles_found,
        'high_risk_diseases': high_risk_diseases,
        'moderate_risk_diseases': moderate_risk_diseases,
        'low_risk_diseases': low_risk_diseases
    }

if __name__ == "__main__":
    analyze_disease_risk() 