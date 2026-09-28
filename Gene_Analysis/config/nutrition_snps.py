"""Nutrition and vitamin metabolism SNP database."""

from typing import Dict

# Nutrition and vitamin metabolism markers
NUTRITION_SNPS = {
    # Vitamin D metabolism
    'rs2282679': {
        'trait': 'Vitamin D levels',
        'gene': 'GC (Vitamin D Binding Protein)',
        'risk_allele': 'C',
        'description': 'Associated with lower vitamin D levels',
        'recommendation': 'May need higher vitamin D supplementation'
    },
    'rs10741657': {
        'trait': 'Vitamin D synthesis',
        'gene': 'CYP2R1',
        'risk_allele': 'A',
        'description': 'Reduced vitamin D synthesis',
        'recommendation': 'Consider vitamin D supplementation'
    },
    'rs12785878': {
        'trait': 'Vitamin D levels',
        'gene': 'NADSYN1',
        'risk_allele': 'T',
        'description': 'Associated with lower vitamin D levels',
        'recommendation': 'Monitor vitamin D levels'
    },

    # Folate/Vitamin B9 metabolism (MTHFR)
    'rs1801133': {
        'trait': 'Folate metabolism',
        'gene': 'MTHFR',
        'variant': 'C677T',
        'risk_allele': 'T',
        'description': 'Reduced MTHFR activity - impaired folate conversion',
        'recommendation': 'Consider methylfolate supplementation, avoid folic acid',
        'clinical_impact': 'Increased homocysteine levels'
    },
    'rs1801131': {
        'trait': 'Folate metabolism',
        'gene': 'MTHFR',
        'variant': 'A1298C',
        'risk_allele': 'C',
        'description': 'Mildly reduced MTHFR activity',
        'recommendation': 'Monitor folate levels'
    },

    # Vitamin B12 metabolism
    'rs602662': {
        'trait': 'Vitamin B12 levels',
        'gene': 'FUT2',
        'risk_allele': 'G',
        'description': 'Associated with lower B12 levels',
        'recommendation': 'Monitor B12 levels, may need supplementation'
    },
    'rs1801222': {
        'trait': 'Vitamin B12 absorption',
        'gene': 'TCN2',
        'risk_allele': 'G',
        'description': 'Reduced B12 transport',
        'recommendation': 'May need higher B12 intake'
    },

    # Omega-3 metabolism
    'rs174537': {
        'trait': 'Omega-3 conversion',
        'gene': 'FADS1',
        'risk_allele': 'G',
        'description': 'Reduced conversion of ALA to EPA/DHA',
        'recommendation': 'Consider direct EPA/DHA supplementation'
    },
    'rs174575': {
        'trait': 'Omega-3 metabolism',
        'gene': 'FADS2',
        'risk_allele': 'C',
        'description': 'Reduced omega-3 fatty acid metabolism',
        'recommendation': 'May benefit from pre-formed omega-3s'
    },

    # Caffeine metabolism (already in fitness, but nutrition-related)
    'rs762551': {
        'trait': 'Caffeine metabolism',
        'gene': 'CYP1A2',
        'risk_allele': 'A',
        'description': 'Slow caffeine metabolism',
        'recommendation': 'Limit caffeine intake, longer effects'
    },

    # Lactose intolerance (already in disease/ancestry)
    'rs4988235': {
        'trait': 'Lactose tolerance',
        'gene': 'LCT',
        'risk_allele': 'C',
        'description': 'Lactose intolerance',
        'recommendation': 'Avoid or limit dairy, consider lactase supplements'
    },

    # Iron metabolism
    'rs1800562': {
        'trait': 'Iron overload',
        'gene': 'HFE',
        'risk_allele': 'G',
        'description': 'Hemochromatosis risk - iron overload',
        'recommendation': 'Limit iron intake, avoid iron supplements'
    },

    # Salt sensitivity
    'rs699': {
        'trait': 'Salt sensitivity',
        'gene': 'AGT',
        'risk_allele': 'A',
        'description': 'Increased salt sensitivity',
        'recommendation': 'Limit sodium intake for blood pressure control'
    },

    # Choline metabolism
    'rs7946': {
        'trait': 'Choline requirements',
        'gene': 'PEMT',
        'risk_allele': 'G',
        'description': 'Increased choline requirements',
        'recommendation': 'Ensure adequate choline intake (eggs, liver, soy)'
    },
}

# Vitamin recommendations based on genotypes
VITAMIN_RECOMMENDATIONS = {
    'Vitamin D': {
        'snps': ['rs2282679', 'rs10741657', 'rs12785878'],
        'general': 'Consider 2000-4000 IU daily if variants present'
    },
    'Folate': {
        'snps': ['rs1801133', 'rs1801131'],
        'general': 'Use methylfolate instead of folic acid if MTHFR variants present'
    },
    'B12': {
        'snps': ['rs602662', 'rs1801222'],
        'general': 'Monitor levels, may need supplementation'
    },
}
