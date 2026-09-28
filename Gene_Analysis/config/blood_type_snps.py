"""Blood type genetic markers (ABO and RHD)."""

from typing import Dict, Tuple

# ABO blood group system
ABO_SNPS = {
    'rs8176719': {
        'system': 'ABO',
        'allele': 'O',
        'genotype_interpretation': {
            '--': 'Blood type O',
            'A-': 'Blood type A',
            'B-': 'Blood type B',
            'AB': 'Blood type AB'
        },
        'description': 'ABO blood group - deletion variant for O allele'
    },
    'rs8176746': {
        'system': 'ABO',
        'allele': 'B',
        'description': 'ABO blood group - B allele marker'
    },
    'rs8176747': {
        'system': 'ABO',
        'allele': 'A',
        'description': 'ABO blood group - A allele marker'
    },
}

# RHD (Rh factor) system
RHD_SNPS = {
    'rs590787': {
        'system': 'RHD',
        'allele': 'Rh+',
        'description': 'RHD positive (Rh+)',
        'genotype_interpretation': {
            'CC': 'Rh positive',
            'CT': 'Rh positive',
            'TT': 'Rh negative'
        }
    },
    'rs590988': {
        'system': 'RHD',
        'allele': 'Rh-',
        'description': 'RHD negative (Rh-) marker'
    },
}

# Blood type combinations
BLOOD_TYPE_COMBINATIONS = {
    'O+': {'ABO': 'O', 'RHD': '+'},
    'O-': {'ABO': 'O', 'RHD': '-'},
    'A+': {'ABO': 'A', 'RHD': '+'},
    'A-': {'ABO': 'A', 'RHD': '-'},
    'B+': {'ABO': 'B', 'RHD': '+'},
    'B-': {'ABO': 'B', 'RHD': '-'},
    'AB+': {'ABO': 'AB', 'RHD': '+'},
    'AB-': {'ABO': 'AB', 'RHD': '-'},
}

# Blood type compatibility (for reference)
BLOOD_TYPE_COMPATIBILITY = {
    'O-': {'can_receive': ['O-'], 'can_donate_to': ['O+', 'O-', 'A+', 'A-', 'B+', 'B-', 'AB+', 'AB-']},
    'O+': {'can_receive': ['O+', 'O-'], 'can_donate_to': ['O+', 'A+', 'B+', 'AB+']},
    'A-': {'can_receive': ['O-', 'A-'], 'can_donate_to': ['A+', 'A-', 'AB+', 'AB-']},
    'A+': {'can_receive': ['O+', 'O-', 'A+', 'A-'], 'can_donate_to': ['A+', 'AB+']},
    'B-': {'can_receive': ['O-', 'B-'], 'can_donate_to': ['B+', 'B-', 'AB+', 'AB-']},
    'B+': {'can_receive': ['O+', 'O-', 'B+', 'B-'], 'can_donate_to': ['B+', 'AB+']},
    'AB-': {'can_receive': ['O-', 'A-', 'B-', 'AB-'], 'can_donate_to': ['AB+', 'AB-']},
    'AB+': {'can_receive': ['O+', 'O-', 'A+', 'A-', 'B+', 'B-', 'AB+', 'AB-'], 'can_donate_to': ['AB+']},
}
