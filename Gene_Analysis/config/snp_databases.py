"""Centralized SNP database with all genetic markers organized by category."""

from typing import Dict, List, Tuple, Any


# ============================================================================
# DISEASE RISK SNPS
# ============================================================================

DISEASE_SNPS = {
    # Global/common disease markers
    'rs429358': {
        'disease': "Alzheimer's (late onset)",
        'risk': ['C'],
        'description': "APOE ε4 allele increases risk of late-onset Alzheimer's.",
        'population': 'Global',
        'gene': 'APOE'
    },
    'rs2187668': {
        'disease': "Celiac disease",
        'risk': ['T'],
        'description': "Associated with HLA-DQ2, key in immune response to gluten.",
        'population': 'Global',
        'gene': 'HLA-DQ2'
    },
    'rs4988235': {
        'disease': "Lactose intolerance",
        'risk': ['C'],
        'description': "CC genotype likely causes lactose intolerance.",
        'population': 'Global',
        'gene': 'LCT'
    },
    'rs1800562': {
        'disease': "Hemochromatosis",
        'risk': ['G'],
        'description': "Mutation in HFE gene leads to iron overload.",
        'population': 'Global',
        'gene': 'HFE'
    },
    'rs6025': {
        'disease': "Factor V Leiden (thrombophilia)",
        'risk': ['A'],
        'description': "Increased risk of blood clots.",
        'population': 'Global',
        'gene': 'F5'
    },
    'rs7903146': {
        'disease': "Type 2 Diabetes",
        'risk': ['T'],
        'description': "TCF7L2 gene variant raises diabetes risk.",
        'population': 'Global',
        'gene': 'TCF7L2'
    },
    'rs1333049': {
        'disease': "Coronary artery disease",
        'risk': ['C'],
        'description': "Strong association with heart disease.",
        'population': 'Global',
        'gene': '9p21'
    },
    'rs10490924': {
        'disease': "Macular degeneration",
        'risk': ['T'],
        'description': "Increases risk of age-related vision loss.",
        'population': 'Global',
        'gene': 'ARMS2'
    },
    'rs2066847': {
        'disease': "Crohn's disease",
        'risk': ['C'],
        'description': "Mutation in NOD2 gene linked to Crohn's.",
        'population': 'Global',
        'gene': 'NOD2'
    },
    'rs2476601': {
        'disease': "Rheumatoid arthritis",
        'risk': ['A'],
        'description': "PTPN22 gene variant increases autoimmune risk.",
        'population': 'Global',
        'gene': 'PTPN22'
    },
    'rs3135388': {
        'disease': "Multiple sclerosis",
        'risk': ['T'],
        'description': "HLA-DRB1*15:01 allele linked to MS.",
        'population': 'Global',
        'gene': 'HLA-DRB1'
    },
    'rs9939609': {
        'disease': "Obesity",
        'risk': ['A'],
        'description': "FTO gene variant associated with higher BMI.",
        'population': 'Global',
        'gene': 'FTO'
    },
    'rs356219': {
        'disease': "Parkinson's disease",
        'risk': ['G'],
        'description': "SNCA gene variant increases risk.",
        'population': 'Global',
        'gene': 'SNCA'
    },
    'rs10484554': {
        'disease': "Psoriasis",
        'risk': ['T'],
        'description': "Linked to immune skin response.",
        'population': 'Global',
        'gene': 'HLA-C'
    },
    'rs10993994': {
        'disease': "Prostate cancer",
        'risk': ['T'],
        'description': "Risk allele in MSMB gene region.",
        'population': 'Global',
        'gene': 'MSMB'
    },
    'rs1799950': {
        'disease': "Breast cancer (BRCA1 proxy)",
        'risk': ['G'],
        'description': "Rare variant possibly linked to BRCA1.",
        'population': 'Global',
        'gene': 'BRCA1'
    },
    'rs6265': {
        'disease': "Depression / neuroticism",
        'risk': ['C'],
        'description': "BDNF gene variant may affect mood regulation.",
        'population': 'Global',
        'gene': 'BDNF'
    },
    'rs12134493': {
        'disease': "Migraine",
        'risk': ['A'],
        'description': "Variant in CACNA1A gene affects migraine susceptibility.",
        'population': 'Global',
        'gene': 'CACNA1A'
    },
    'rs7216389': {
        'disease': "Asthma (childhood)",
        'risk': ['T'],
        'description': "Variant on 17q21 linked to early asthma.",
        'population': 'Global',
        'gene': '17q21'
    },
    'rs7574865': {
        'disease': "Lupus (SLE)",
        'risk': ['T'],
        'description': "STAT4 gene variant common in autoimmunity.",
        'population': 'Global',
        'gene': 'STAT4'
    },
    # East Asian-specific markers
    'rs671': {
        'disease': "Alcohol flush reaction (ALDH2 deficiency)",
        'risk': ['A'],
        'description': "ALDH2*2 allele causes alcohol intolerance, common in East Asians.",
        'population': 'East Asian',
        'gene': 'ALDH2'
    },
    'rs1229984': {
        'disease': "Alcohol metabolism (ADH1B)",
        'risk': ['A'],
        'description': "ADH1B*2 allele increases alcohol metabolism, common in East Asians.",
        'population': 'East Asian',
        'gene': 'ADH1B'
    },
    'rs2075650': {
        'disease': "Alzheimer's disease (APOE region, East Asian)",
        'risk': ['G'],
        'description': "Associated with Alzheimer's in East Asians.",
        'population': 'East Asian',
        'gene': 'APOE'
    },
    'rs1801133': {
        'disease': "Homocysteine metabolism (MTHFR)",
        'risk': ['T'],
        'description': "MTHFR C677T variant, higher risk of hyperhomocysteinemia, common in East Asians.",
        'population': 'East Asian',
        'gene': 'MTHFR'
    },
    'rs2231142': {
        'disease': "Gout (ABCG2)",
        'risk': ['T'],
        'description': "ABCG2 Q141K variant, high gout risk in East Asians.",
        'population': 'East Asian',
        'gene': 'ABCG2'
    },
    'rs2285666': {
        'disease': "ACE2 expression (COVID-19 susceptibility)",
        'risk': ['A'],
        'description': "Variant may affect ACE2 expression, studied in East Asians.",
        'population': 'East Asian',
        'gene': 'ACE2'
    },
    'rs11200638': {
        'disease': "Age-related macular degeneration (HTRA1)",
        'risk': ['A'],
        'description': "HTRA1 risk allele, high prevalence in East Asians.",
        'population': 'East Asian',
        'gene': 'HTRA1'
    },
    'rs1800414': {
        'disease': "Skin pigmentation (OCA2)",
        'risk': ['G'],
        'description': "OCA2 variant, common in East Asians.",
        'population': 'East Asian',
        'gene': 'OCA2'
    },
    'rs1042522': {
        'disease': "Cancer risk (TP53)",
        'risk': ['C'],
        'description': "TP53 Arg72Pro, cancer risk variant, higher in East Asians.",
        'population': 'East Asian',
        'gene': 'TP53'
    },
    'rs11655237': {
        'disease': "Gastric cancer (LINC00673)",
        'risk': ['A'],
        'description': "LINC00673 variant, gastric cancer risk in East Asians.",
        'population': 'East Asian',
        'gene': 'LINC00673'
    },
    'rs2736100': {
        'disease': "Lung cancer (TERT)",
        'risk': ['A'],
        'description': "TERT variant, lung cancer risk in East Asians.",
        'population': 'East Asian',
        'gene': 'TERT'
    },
    'rs1801274': {
        'disease': "Autoimmune disease (FCGR2A)",
        'risk': ['A'],
        'description': "FCGR2A variant, SLE/autoimmunity risk, higher in East Asians.",
        'population': 'East Asian',
        'gene': 'FCGR2A'
    },
}

# Disease statistics (population frequencies, odds ratios, etc.)
DISEASE_STATS = {
    'rs429358': {
        'gene': 'APOE',
        'population_frequency': 0.15,
        'odds_ratio': 3.5,
        'confidence_interval': '2.8-4.4'
    },
    'rs2187668': {
        'gene': 'HLA-DQ2',
        'population_frequency': 0.08,
        'odds_ratio': 2.8,
        'confidence_interval': '2.1-3.7'
    },
    'rs4988235': {
        'gene': 'LCT',
        'population_frequency': 0.65,
        'odds_ratio': 4.2,
        'confidence_interval': '3.5-5.1'
    },
    'rs1800562': {
        'gene': 'HFE',
        'population_frequency': 0.06,
        'odds_ratio': 8.5,
        'confidence_interval': '6.2-11.8'
    },
    'rs6025': {
        'gene': 'F5',
        'population_frequency': 0.05,
        'odds_ratio': 3.8,
        'confidence_interval': '2.9-5.0'
    },
    'rs7903146': {
        'gene': 'TCF7L2',
        'population_frequency': 0.30,
        'odds_ratio': 1.4,
        'confidence_interval': '1.2-1.6'
    },
    'rs1333049': {
        'gene': '9p21',
        'population_frequency': 0.45,
        'odds_ratio': 1.3,
        'confidence_interval': '1.1-1.5'
    },
    'rs10490924': {
        'gene': 'ARMS2',
        'population_frequency': 0.25,
        'odds_ratio': 2.1,
        'confidence_interval': '1.7-2.6'
    },
    'rs2066847': {
        'gene': 'NOD2',
        'population_frequency': 0.10,
        'odds_ratio': 2.5,
        'confidence_interval': '1.9-3.3'
    },
    'rs2476601': {
        'gene': 'PTPN22',
        'population_frequency': 0.12,
        'odds_ratio': 1.8,
        'confidence_interval': '1.4-2.3'
    },
    'rs3135388': {
        'gene': 'HLA-DRB1',
        'population_frequency': 0.15,
        'odds_ratio': 2.2,
        'confidence_interval': '1.7-2.8'
    },
    'rs9939609': {
        'gene': 'FTO',
        'population_frequency': 0.40,
        'odds_ratio': 1.2,
        'confidence_interval': '1.1-1.3'
    },
    'rs356219': {
        'gene': 'SNCA',
        'population_frequency': 0.20,
        'odds_ratio': 1.4,
        'confidence_interval': '1.1-1.7'
    },
    'rs10484554': {
        'gene': 'HLA-C',
        'population_frequency': 0.18,
        'odds_ratio': 1.6,
        'confidence_interval': '1.3-2.0'
    },
    'rs10993994': {
        'gene': 'MSMB',
        'population_frequency': 0.35,
        'odds_ratio': 1.3,
        'confidence_interval': '1.1-1.5'
    },
    'rs1799950': {
        'gene': 'BRCA1',
        'population_frequency': 0.02,
        'odds_ratio': 2.8,
        'confidence_interval': '1.8-4.3'
    },
    'rs6265': {
        'gene': 'BDNF',
        'population_frequency': 0.25,
        'odds_ratio': 1.2,
        'confidence_interval': '1.0-1.4'
    },
    'rs12134493': {
        'gene': 'CACNA1A',
        'population_frequency': 0.30,
        'odds_ratio': 1.3,
        'confidence_interval': '1.1-1.5'
    },
    'rs7216389': {
        'gene': '17q21',
        'population_frequency': 0.22,
        'odds_ratio': 1.4,
        'confidence_interval': '1.2-1.7'
    },
    'rs7574865': {
        'gene': 'STAT4',
        'population_frequency': 0.15,
        'odds_ratio': 1.6,
        'confidence_interval': '1.3-2.0'
    },
}

# Polygenic Risk Score (PRS) models
DISEASE_PRS_MODELS = {
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


# ============================================================================
# LONGEVITY MARKERS
# ============================================================================

LONGEVITY_MARKERS = {
    # Major markers
    'rs2736100': {
        'risk': ['A'],
        'weight': 0.15,
        'description': 'TERT – Telomere maintenance'
    },
    'rs7726159': {
        'risk': ['A'],
        'weight': 0.12,
        'description': 'TERT – Telomere length'
    },
    'rs1317082': {
        'risk': ['A'],
        'weight': 0.10,
        'description': 'TERC – Telomere regulation'
    },
    'rs1801133': {
        'risk': ['A'],
        'weight': 0.19,
        'description': 'MTHFR – Methylation cycle'
    },
    'rs1042522': {
        'risk': ['C'],
        'weight': 0.10,
        'description': 'TP53 – DNA repair'
    },
    # Additional research markers
    'rs10936599': {
        'risk': ['C'],
        'weight': 0.05,
        'description': 'TERC – Telomere length'
    },
    'rs755017': {
        'risk': ['A'],
        'weight': 0.04,
        'description': 'FOXO3 – Longevity variant'
    },
    'rs2802292': {
        'risk': ['G'],
        'weight': 0.03,
        'description': 'FOXO3 – Longevity variant'
    },
    'rs2157719': {
        'risk': ['C'],
        'weight': 0.03,
        'description': 'CDKN2B – Cellular senescence'
    },
    'rs11125529': {
        'risk': ['T'],
        'weight': 0.02,
        'description': 'ACYP2 – Telomere length'
    },
}

# Total markers for normalization (used in calculations)
LONGEVITY_TOTAL_MARKERS = 25


# ============================================================================
# FITNESS & ATHLETIC PERFORMANCE SNPS
# ============================================================================

FITNESS_SNPS = {
    'Muscle Fiber Type': {
        'rs1815739': {'risk': ['C'], 'description': 'ACTN3 - Fast twitch muscle fibers', 'weight': 0.30},
        'rs540874': {'risk': ['A'], 'description': 'ACE - Endurance performance', 'weight': 0.25},
        'rs4343': {'risk': ['G'], 'description': 'ACE - Power performance', 'weight': 0.20},
        'rs699': {'risk': ['A'], 'description': 'AGT - Blood pressure regulation', 'weight': 0.15},
        'rs5186': {'risk': ['C'], 'description': 'AGTR1 - Muscle hypertrophy', 'weight': 0.10},
    },
    'Endurance Capacity': {
        'rs1815739': {'risk': ['T'], 'description': 'ACTN3 - Endurance advantage', 'weight': 0.35},
        'rs540874': {'risk': ['A'], 'description': 'ACE - Endurance performance', 'weight': 0.30},
        'rs4343': {'risk': ['A'], 'description': 'ACE - Endurance variant', 'weight': 0.20},
        'rs2070744': {'risk': ['C'], 'description': 'NOS3 - Nitric oxide production', 'weight': 0.15},
    },
    'Power & Strength': {
        'rs1815739': {'risk': ['C'], 'description': 'ACTN3 - Power performance', 'weight': 0.40},
        'rs4343': {'risk': ['G'], 'description': 'ACE - Power variant', 'weight': 0.30},
        'rs5186': {'risk': ['C'], 'description': 'AGTR1 - Muscle strength', 'weight': 0.20},
        'rs699': {'risk': ['G'], 'description': 'AGT - Strength performance', 'weight': 0.10},
    },
    'Injury Risk': {
        'rs1800012': {'risk': ['A'], 'description': 'COL1A1 - ACL injury risk', 'weight': 0.35},
        'rs12722': {'risk': ['T'], 'description': 'COL5A1 - Tendon injury risk', 'weight': 0.30},
        'rs13946': {'risk': ['C'], 'description': 'COL12A1 - Ligament injury', 'weight': 0.20},
        'rs1800013': {'risk': ['G'], 'description': 'COL1A1 - Stress fracture risk', 'weight': 0.15},
    },
    'Recovery Rate': {
        'rs1800629': {'risk': ['A'], 'description': 'TNF - Inflammation response', 'weight': 0.30},
        'rs1801133': {'risk': ['T'], 'description': 'MTHFR - Recovery factors', 'weight': 0.25},
        'rs1801131': {'risk': ['C'], 'description': 'MTHFR - Recovery metabolism', 'weight': 0.20},
        'rs2070744': {'risk': ['T'], 'description': 'NOS3 - Recovery efficiency', 'weight': 0.15},
        'rs1799983': {'risk': ['T'], 'description': 'NOS3 - Recovery optimization', 'weight': 0.10},
    },
    'VO2 Max Potential': {
        'rs1815739': {'risk': ['T'], 'description': 'ACTN3 - Aerobic capacity', 'weight': 0.35},
        'rs540874': {'risk': ['A'], 'description': 'ACE - VO2 max', 'weight': 0.30},
        'rs4343': {'risk': ['A'], 'description': 'ACE - Endurance capacity', 'weight': 0.20},
        'rs2070744': {'risk': ['C'], 'description': 'NOS3 - Oxygen utilization', 'weight': 0.15},
    },
    'Muscle Hypertrophy': {
        'rs1815739': {'risk': ['C'], 'description': 'ACTN3 - Muscle growth', 'weight': 0.30},
        'rs5186': {'risk': ['C'], 'description': 'AGTR1 - Muscle hypertrophy', 'weight': 0.25},
        'rs4343': {'risk': ['G'], 'description': 'ACE - Muscle development', 'weight': 0.20},
        'rs699': {'risk': ['G'], 'description': 'AGT - Muscle building', 'weight': 0.15},
        'rs1800629': {'risk': ['G'], 'description': 'TNF - Muscle response', 'weight': 0.10},
    },
    'Lactate Threshold': {
        'rs1815739': {'risk': ['T'], 'description': 'ACTN3 - Lactate clearance', 'weight': 0.35},
        'rs540874': {'risk': ['A'], 'description': 'ACE - Lactate threshold', 'weight': 0.30},
        'rs4343': {'risk': ['A'], 'description': 'ACE - Endurance efficiency', 'weight': 0.20},
        'rs2070744': {'risk': ['C'], 'description': 'NOS3 - Metabolic efficiency', 'weight': 0.15},
    },
    'Caffeine Response': {
        'rs762551': {'risk': ['A'], 'description': 'CYP1A2 - Caffeine metabolism', 'weight': 0.40},
        'rs5751876': {'risk': ['T'], 'description': 'ADORA2A - Caffeine sensitivity', 'weight': 0.30},
        'rs5751879': {'risk': ['C'], 'description': 'ADORA2A - Performance response', 'weight': 0.20},
        'rs2470893': {'risk': ['T'], 'description': 'CYP1A2 - Caffeine effects', 'weight': 0.10},
    },
    'Heat Tolerance': {
        'rs1800629': {'risk': ['A'], 'description': 'TNF - Heat stress response', 'weight': 0.35},
        'rs1801133': {'risk': ['T'], 'description': 'MTHFR - Heat adaptation', 'weight': 0.25},
        'rs2070744': {'risk': ['C'], 'description': 'NOS3 - Thermoregulation', 'weight': 0.20},
        'rs1799983': {'risk': ['T'], 'description': 'NOS3 - Heat tolerance', 'weight': 0.20},
    }
}


# ============================================================================
# ANCESTRY MARKERS
# ============================================================================

ANCESTRY_CATEGORIES = {
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
    'European Ancestry Markers': {
        'rs6058017': 'HERC2 - European light pigmentation',
    }
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_all_snps_by_category() -> Dict[str, List[str]]:
    """Get all SNPs organized by category."""
    categories = {
        'disease': list(DISEASE_SNPS.keys()),
        'longevity': list(LONGEVITY_MARKERS.keys()),
        'fitness': [],
        'ancestry': [],
    }

    # Collect fitness SNPs
    for trait_snps in FITNESS_SNPS.values():
        categories['fitness'].extend(trait_snps.keys())
    categories['fitness'] = list(set(categories['fitness']))  # Remove duplicates

    # Collect ancestry SNPs
    for category_snps in ANCESTRY_CATEGORIES.values():
        categories['ancestry'].extend(category_snps.keys())
    categories['ancestry'] = list(set(categories['ancestry']))  # Remove duplicates

    return categories
