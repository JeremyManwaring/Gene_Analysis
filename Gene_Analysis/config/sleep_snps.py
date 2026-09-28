"""Sleep and chronotype genetic markers."""

from typing import Dict

# Sleep and chronotype markers
SLEEP_SNPS = {
    # Circadian rhythm genes
    'rs1801260': {
        'trait': 'Chronotype (morning/evening preference)',
        'gene': 'CLOCK',
        'risk_allele': 'C',
        'description': 'Associated with evening preference (night owl)',
        'interpretation': 'May prefer later sleep/wake times'
    },
    'rs10462021': {
        'trait': 'Chronotype',
        'gene': 'PER2',
        'risk_allele': 'G',
        'description': 'Associated with morning preference',
        'interpretation': 'May prefer earlier sleep/wake times'
    },
    'rs228697': {
        'trait': 'Chronotype',
        'gene': 'PER3',
        'risk_allele': 'C',
        'description': 'Associated with evening preference',
        'interpretation': 'Evening chronotype'
    },

    # Sleep duration
    'rs11046212': {
        'trait': 'Sleep duration',
        'gene': 'PAX8',
        'risk_allele': 'T',
        'description': 'Associated with shorter sleep duration',
        'interpretation': 'May naturally need less sleep'
    },

    # Sleep quality
    'rs73598374': {
        'trait': 'Sleep quality',
        'gene': 'DEC2',
        'risk_allele': 'C',
        'description': 'Natural short sleeper variant',
        'interpretation': 'May function well on less sleep'
    },

    # Insomnia risk
    'rs113851554': {
        'trait': 'Insomnia risk',
        'gene': 'MEIS1',
        'risk_allele': 'C',
        'description': 'Increased insomnia risk',
        'interpretation': 'Higher genetic predisposition to insomnia'
    },

    # Restless leg syndrome
    'rs3923809': {
        'trait': 'Restless leg syndrome',
        'gene': 'BTBD9',
        'risk_allele': 'A',
        'description': 'Increased restless leg syndrome risk',
        'interpretation': 'May experience RLS symptoms'
    },

    # Narcolepsy (HLA association)
    'rs1154155': {
        'trait': 'Narcolepsy risk',
        'gene': 'HLA-DQB1',
        'risk_allele': 'C',
        'description': 'Narcolepsy risk marker',
        'interpretation': 'Genetic risk factor for narcolepsy'
    },

    # Caffeine sensitivity (affects sleep)
    'rs5751876': {
        'trait': 'Caffeine sensitivity',
        'gene': 'ADORA2A',
        'risk_allele': 'T',
        'description': 'Increased caffeine sensitivity - affects sleep more',
        'interpretation': 'Caffeine may significantly disrupt sleep'
    },
}

# Chronotype interpretation
CHRONOTYPE_INTERPRETATION = {
    'Morning': 'Early riser, peak performance in morning',
    'Evening': 'Night owl, peak performance in evening',
    'Intermediate': 'Flexible sleep schedule'
}
