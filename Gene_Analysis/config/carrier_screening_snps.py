"""Carrier screening genetic markers for inherited conditions."""

from typing import Dict

# Carrier screening markers
CARRIER_SCREENING_SNPS = {
    # Cystic Fibrosis (CFTR gene)
    'rs113993960': {
        'condition': 'Cystic Fibrosis',
        'gene': 'CFTR',
        'variant': 'F508del',
        'risk_allele': 'del',
        'description': 'Most common CFTR mutation',
        'carrier_status': 'Carrier if heterozygous',
        'prevalence': '1 in 25-30 Caucasians are carriers'
    },
    'rs75527207': {
        'condition': 'Cystic Fibrosis',
        'gene': 'CFTR',
        'variant': 'G542X',
        'risk_allele': 'G',
        'description': 'CFTR nonsense mutation',
        'carrier_status': 'Carrier if heterozygous'
    },
    'rs121908755': {
        'condition': 'Cystic Fibrosis',
        'gene': 'CFTR',
        'variant': 'G551D',
        'risk_allele': 'G',
        'description': 'CFTR gating mutation',
        'carrier_status': 'Carrier if heterozygous'
    },

    # Sickle Cell Anemia (HBB gene)
    'rs334': {
        'condition': 'Sickle Cell Anemia',
        'gene': 'HBB',
        'variant': 'HbS',
        'risk_allele': 'T',
        'description': 'Sickle cell mutation',
        'carrier_status': 'Carrier (sickle cell trait) if heterozygous',
        'prevalence': '1 in 13 African Americans are carriers'
    },

    # Beta Thalassemia (HBB gene)
    'rs11549407': {
        'condition': 'Beta Thalassemia',
        'gene': 'HBB',
        'variant': 'Beta-thalassemia',
        'risk_allele': 'C',
        'description': 'Beta thalassemia mutation',
        'carrier_status': 'Carrier if heterozygous'
    },

    # Tay-Sachs Disease (HEXA gene)
    'rs121907957': {
        'condition': 'Tay-Sachs Disease',
        'gene': 'HEXA',
        'variant': '1278+1G>C',
        'risk_allele': 'G',
        'description': 'Tay-Sachs mutation',
        'carrier_status': 'Carrier if heterozygous',
        'prevalence': '1 in 27 Ashkenazi Jews are carriers'
    },

    # Spinal Muscular Atrophy (SMN1 gene)
    'rs121908192': {
        'condition': 'Spinal Muscular Atrophy',
        'gene': 'SMN1',
        'variant': 'SMN1 deletion',
        'risk_allele': 'del',
        'description': 'SMN1 gene deletion',
        'carrier_status': 'Carrier if heterozygous',
        'prevalence': '1 in 40-50 people are carriers'
    },

    # Fragile X Syndrome (FMR1 gene - CGG repeat, not SNP)
    # Note: This requires special testing, not standard SNP array

    # Gaucher Disease (GBA gene)
    'rs421016': {
        'condition': 'Gaucher Disease',
        'gene': 'GBA',
        'variant': 'N370S',
        'risk_allele': 'A',
        'description': 'Gaucher disease mutation',
        'carrier_status': 'Carrier if heterozygous',
        'prevalence': '1 in 15 Ashkenazi Jews are carriers'
    },

    # Canavan Disease (ASPA gene)
    'rs121907957': {
        'condition': 'Canavan Disease',
        'gene': 'ASPA',
        'variant': 'E285A',
        'risk_allele': 'A',
        'description': 'Canavan disease mutation',
        'carrier_status': 'Carrier if heterozygous'
    },

    # Hemochromatosis (already in disease, but also carrier screening)
    'rs1800562': {
        'condition': 'Hemochromatosis',
        'gene': 'HFE',
        'variant': 'C282Y',
        'risk_allele': 'G',
        'description': 'Hemochromatosis mutation',
        'carrier_status': 'Carrier if heterozygous',
        'prevalence': '1 in 10 people of Northern European descent are carriers'
    },

    # Alpha-1 Antitrypsin Deficiency (SERPINA1 gene)
    'rs28929474': {
        'condition': 'Alpha-1 Antitrypsin Deficiency',
        'gene': 'SERPINA1',
        'variant': 'Pi*Z',
        'risk_allele': 'G',
        'description': 'Alpha-1 antitrypsin deficiency',
        'carrier_status': 'Carrier if heterozygous'
    },
}

# Population-specific carrier frequencies
POPULATION_CARRIER_FREQUENCIES = {
    'Cystic Fibrosis': {
        'Caucasian': '1 in 25-30',
        'Hispanic': '1 in 46',
        'African American': '1 in 65',
        'Asian': '1 in 90'
    },
    'Sickle Cell Anemia': {
        'African American': '1 in 13',
        'Hispanic': '1 in 140',
        'Caucasian': '1 in 1000'
    },
    'Tay-Sachs Disease': {
        'Ashkenazi Jewish': '1 in 27',
        'French Canadian': '1 in 30',
        'Cajun': '1 in 50',
        'General population': '1 in 250'
    },
}
