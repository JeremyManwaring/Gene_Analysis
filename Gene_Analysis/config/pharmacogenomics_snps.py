"""Pharmacogenomics SNP database - drug metabolism and response markers."""

from typing import Dict, Tuple

# Pharmacogenomics markers for drug metabolism
PHARMACOGENOMICS_SNPS = {
    # CYP2D6 - Metabolizes ~25% of drugs (antidepressants, antipsychotics, beta-blockers)
    'rs1065852': {
        'gene': 'CYP2D6',
        'drugs': ['Codeine', 'Tramadol', 'Amitriptyline', 'Metoprolol'],
        'variant': 'CYP2D6*10 (reduced function)',
        'risk_allele': 'C',
        'description': 'Reduced CYP2D6 enzyme activity - may need lower doses',
        'phenotype': 'Intermediate metabolizer'
    },
    'rs16947': {
        'gene': 'CYP2D6',
        'drugs': ['Codeine', 'Tramadol', 'Amitriptyline'],
        'variant': 'CYP2D6*2 (normal function)',
        'risk_allele': 'G',
        'description': 'Normal CYP2D6 enzyme activity',
        'phenotype': 'Normal metabolizer'
    },
    'rs35742686': {
        'gene': 'CYP2D6',
        'drugs': ['Codeine', 'Tramadol'],
        'variant': 'CYP2D6*4 (no function)',
        'risk_allele': 'G',
        'description': 'No CYP2D6 enzyme activity - poor metabolizer',
        'phenotype': 'Poor metabolizer'
    },

    # CYP2C19 - Metabolizes antidepressants, antiplatelets, proton pump inhibitors
    'rs4244285': {
        'gene': 'CYP2C19',
        'drugs': ['Clopidogrel', 'Omeprazole', 'Sertraline', 'Citalopram'],
        'variant': 'CYP2C19*2 (no function)',
        'risk_allele': 'A',
        'description': 'No CYP2C19 enzyme activity - poor metabolizer',
        'phenotype': 'Poor metabolizer',
        'clinical_impact': 'Clopidogrel may be less effective'
    },
    'rs4986893': {
        'gene': 'CYP2C19',
        'drugs': ['Clopidogrel', 'Omeprazole'],
        'variant': 'CYP2C19*3 (no function)',
        'risk_allele': 'G',
        'description': 'No CYP2C19 enzyme activity',
        'phenotype': 'Poor metabolizer'
    },
    'rs12248560': {
        'gene': 'CYP2C19',
        'drugs': ['Omeprazole', 'Pantoprazole'],
        'variant': 'CYP2C19*17 (ultra-rapid)',
        'risk_allele': 'C',
        'description': 'Increased CYP2C19 enzyme activity - ultra-rapid metabolizer',
        'phenotype': 'Ultra-rapid metabolizer'
    },

    # CYP2C9 - Metabolizes warfarin, NSAIDs, phenytoin
    'rs1799853': {
        'gene': 'CYP2C9',
        'drugs': ['Warfarin', 'Ibuprofen', 'Phenytoin'],
        'variant': 'CYP2C9*2 (reduced function)',
        'risk_allele': 'C',
        'description': 'Reduced CYP2C9 activity - may need lower warfarin doses',
        'phenotype': 'Intermediate metabolizer'
    },
    'rs1057910': {
        'gene': 'CYP2C9',
        'drugs': ['Warfarin', 'Ibuprofen'],
        'variant': 'CYP2C9*3 (reduced function)',
        'risk_allele': 'A',
        'description': 'Reduced CYP2C9 activity - significantly lower warfarin doses needed',
        'phenotype': 'Intermediate metabolizer'
    },

    # VKORC1 - Warfarin sensitivity
    'rs9923231': {
        'gene': 'VKORC1',
        'drugs': ['Warfarin'],
        'variant': 'VKORC1 -1639G>A',
        'risk_allele': 'A',
        'description': 'Increased warfarin sensitivity - lower doses needed',
        'phenotype': 'Sensitive'
    },

    # TPMT - Thiopurine metabolism (azathioprine, mercaptopurine)
    'rs1800462': {
        'gene': 'TPMT',
        'drugs': ['Azathioprine', 'Mercaptopurine', 'Thioguanine'],
        'variant': 'TPMT*3C (no function)',
        'risk_allele': 'A',
        'description': 'No TPMT activity - high risk of toxicity, avoid or use very low doses',
        'phenotype': 'Poor metabolizer'
    },
    'rs1142345': {
        'gene': 'TPMT',
        'drugs': ['Azathioprine', 'Mercaptopurine'],
        'variant': 'TPMT*3A (no function)',
        'risk_allele': 'G',
        'description': 'No TPMT activity - high risk of toxicity',
        'phenotype': 'Poor metabolizer'
    },

    # DPYD - Fluoropyrimidine metabolism (5-FU, capecitabine)
    'rs3918290': {
        'gene': 'DPYD',
        'drugs': ['5-Fluorouracil', 'Capecitabine'],
        'variant': 'DPYD*2A (no function)',
        'risk_allele': 'A',
        'description': 'No DPYD activity - severe toxicity risk, avoid or use reduced doses',
        'phenotype': 'Poor metabolizer'
    },
    'rs55886062': {
        'gene': 'DPYD',
        'drugs': ['5-Fluorouracil', 'Capecitabine'],
        'variant': 'DPYD*13 (reduced function)',
        'risk_allele': 'C',
        'description': 'Reduced DPYD activity - may need dose reduction',
        'phenotype': 'Intermediate metabolizer'
    },

    # SLCO1B1 - Statin metabolism
    'rs4149056': {
        'gene': 'SLCO1B1',
        'drugs': ['Simvastatin', 'Atorvastatin', 'Pravastatin'],
        'variant': 'SLCO1B1*5 (reduced function)',
        'risk_allele': 'C',
        'description': 'Reduced statin transport - increased risk of myopathy',
        'phenotype': 'At risk'
    },

    # UGT1A1 - Irinotecan metabolism
    'rs8175347': {
        'gene': 'UGT1A1',
        'drugs': ['Irinotecan'],
        'variant': 'UGT1A1*28 (reduced function)',
        'risk_allele': 'TA7',
        'description': 'Reduced UGT1A1 activity - increased irinotecan toxicity risk',
        'phenotype': 'At risk'
    },
}

# Drug-specific recommendations
DRUG_RECOMMENDATIONS = {
    'Warfarin': {
        'genes': ['CYP2C9', 'VKORC1'],
        'recommendation': 'Dose adjustment based on CYP2C9 and VKORC1 variants'
    },
    'Clopidogrel': {
        'genes': ['CYP2C19'],
        'recommendation': 'Consider alternative antiplatelet if CYP2C19 poor metabolizer'
    },
    'Codeine': {
        'genes': ['CYP2D6'],
        'recommendation': 'Avoid in CYP2D6 poor metabolizers (ineffective) or ultra-rapid (toxicity risk)'
    },
}
