"""Cognitive and intelligence genetic markers."""

from typing import Dict

# Cognitive markers organized by domain
COGNITIVE_SNPS = {
    # Working Memory & Executive Function
    'Working Memory': {
        'rs4680': {
            'gene': 'COMT',
            'variant': 'Val158Met',
            'risk_allele': 'G',  # Val allele
            'description': 'COMT Val158Met - affects dopamine breakdown in prefrontal cortex',
            'interpretation': {
                'GG': 'Val/Val - Better working memory under stress, faster dopamine breakdown',
                'AG': 'Val/Met - Intermediate - Balanced performance',
                'AA': 'Met/Met - Better baseline working memory, slower dopamine breakdown'
            },
            'cognitive_impact': 'Working memory and executive function'
        },
    },

    # Long-term Memory
    'Long-term Memory': {
        'rs17070145': {
            'gene': 'KIBRA',
            'variant': 'KIBRA T/C',
            'risk_allele': 'T',
            'description': 'KIBRA gene variant associated with memory formation and recall',
            'interpretation': {
                'TT': 'Better memory performance',
                'CT': 'Intermediate memory',
                'CC': 'Lower memory performance'
            },
            'cognitive_impact': 'Memory formation and long-term memory'
        },
        'rs6465084': {
            'gene': 'GRM3',
            'variant': 'GRM3 variant',
            'risk_allele': 'A',
            'description': 'GRM3 gene variant associated with memory performance',
            'interpretation': {
                'AA': 'Better memory performance',
                'AG': 'Intermediate memory',
                'GG': 'Lower memory performance'
            },
            'cognitive_impact': 'Memory performance and recall'
        },
        'rs363050': {
            'gene': 'SNAP25',
            'variant': 'SNAP25 variant',
            'risk_allele': 'G',
            'description': 'SNAP25 gene variant associated with memory and learning',
            'interpretation': {
                'GG': 'Better memory and learning',
                'AG': 'Intermediate',
                'AA': 'Lower memory performance'
            },
            'cognitive_impact': 'Memory and synaptic plasticity'
        },
    },

    # Learning Ability
    'Learning Ability': {
        # Neuroplasticity markers
        'rs6265': {
            'gene': 'BDNF',
            'variant': 'Val66Met',
            'risk_allele': 'C',  # Met allele
            'description': 'BDNF Val66Met - brain-derived neurotrophic factor affects learning and memory',
            'interpretation': {
                'GG': 'Val/Val - Better learning and memory, higher BDNF',
                'CG': 'Val/Met - Intermediate learning ability',
                'CC': 'Met/Met - Reduced learning efficiency, lower BDNF secretion'
            },
            'cognitive_impact': 'Learning ability, memory formation, neuroplasticity',
            'category': 'Neuroplasticity'
        },
        'rs7124442': {
            'gene': 'BDNF',
            'variant': 'BDNF promoter',
            'risk_allele': 'A',
            'description': 'BDNF promoter variant - affects BDNF gene expression and production',
            'interpretation': {
                'AA': 'Higher BDNF expression potential',
                'AG': 'Intermediate BDNF expression',
                'GG': 'Standard BDNF expression'
            },
            'cognitive_impact': 'BDNF production, neuroplasticity, learning capacity',
            'category': 'Neuroplasticity'
        },
        # Reward learning markers
        'rs1800497': {
            'gene': 'DRD2',
            'variant': 'Taq1A',
            'risk_allele': 'A',
            'description': 'DRD2 Taq1A variant affects dopamine receptor D2, influences reward learning',
            'interpretation': {
                'AA': 'A1/A1 - Reduced dopamine receptor density, different reward processing',
                'AG': 'A1/A2 - Intermediate reward learning',
                'GG': 'A2/A2 - Higher dopamine receptor density, better reward learning'
            },
            'cognitive_impact': 'Reward processing, learning from feedback, motivation',
            'category': 'Reward Learning'
        },
        'rs6277': {
            'gene': 'DRD2',
            'variant': 'DRD2 C957T',
            'risk_allele': 'C',
            'description': 'DRD2 C957T variant affects dopamine receptor availability and reward learning',
            'interpretation': {
                'CC': 'Higher receptor availability, better reward learning',
                'CT': 'Intermediate reward learning',
                'TT': 'Lower receptor availability'
            },
            'cognitive_impact': 'Reward learning, motivation, reinforcement learning',
            'category': 'Reward Learning'
        },
        'rs1801028': {
            'gene': 'DRD2',
            'variant': 'Ser311Cys',
            'risk_allele': 'G',
            'description': 'DRD2 Ser311Cys variant affects dopamine receptor function',
            'interpretation': {
                'GG': 'Ser/Ser - Standard receptor function',
                'CG': 'Ser/Cys - Altered receptor function',
                'CC': 'Cys/Cys - Different receptor signaling'
            },
            'cognitive_impact': 'Dopamine signaling, learning efficiency',
            'category': 'Reward Learning'
        },
        # Memory consolidation markers
        'rs10835210': {
            'gene': 'CACNA1C',
            'variant': 'CACNA1C variant',
            'risk_allele': 'A',
            'description': 'CACNA1C calcium channel variant affects synaptic plasticity and learning',
            'interpretation': {
                'AA': 'Enhanced synaptic plasticity, better learning',
                'AG': 'Intermediate learning capacity',
                'GG': 'Standard learning capacity'
            },
            'cognitive_impact': 'Synaptic plasticity, memory consolidation, learning',
            'category': 'Memory Consolidation'
        },
        # Attention and focus markers
        'rs25531': {
            'gene': 'SLC6A4',
            'variant': '5-HTTLPR',
            'risk_allele': 'A',
            'description': 'Serotonin transporter variant affects attention, focus, and learning under stress',
            'interpretation': {
                'AA': 'Short/short - Better performance under stress, different attention patterns',
                'AG': 'Short/long - Intermediate attention and focus',
                'GG': 'Long/long - Better baseline attention, may struggle under stress'
            },
            'cognitive_impact': 'Attention, focus, learning under stress, emotional regulation',
            'category': 'Attention/Focus'
        },
        'rs1360780': {
            'gene': 'FKBP5',
            'variant': 'FKBP5 variant',
            'risk_allele': 'T',
            'description': 'FKBP5 variant affects stress response and learning under pressure',
            'interpretation': {
                'TT': 'Enhanced stress response, may affect learning under stress',
                'CT': 'Intermediate stress response',
                'CC': 'Better learning under stress conditions'
            },
            'cognitive_impact': 'Stress response, learning under pressure, cortisol regulation',
            'category': 'Attention/Focus'
        },
    },

    # General Intelligence
    'General Intelligence': {
        'rs324650': {
            'gene': 'CHRM2',
            'variant': 'CHRM2 variant',
            'risk_allele': 'T',
            'description': 'CHRM2 gene variant associated with general cognitive ability',
            'interpretation': {
                'TT': 'Higher general cognitive ability',
                'CT': 'Intermediate',
                'CC': 'Lower general cognitive ability'
            },
            'cognitive_impact': 'General intelligence and cognitive performance'
        },
        'rs7171755': {
            'gene': 'NPTN',
            'variant': 'NPTN variant',
            'risk_allele': 'G',
            'description': 'NPTN gene variant associated with general intelligence',
            'interpretation': {
                'GG': 'Higher general intelligence',
                'AG': 'Intermediate',
                'AA': 'Lower general intelligence'
            },
            'cognitive_impact': 'General intelligence and brain structure'
        },
        'rs174575': {
            'gene': 'FADS2',
            'variant': 'FADS2 variant',
            'risk_allele': 'C',
            'description': 'FADS2 gene variant may interact with nutrition for IQ development',
            'interpretation': {
                'CC': 'May benefit from omega-3 for cognitive development',
                'CT': 'Intermediate',
                'TT': 'Standard cognitive development'
            },
            'cognitive_impact': 'IQ development, may interact with nutrition'
        },
    },

    # Cognitive Decline Risk
    'Cognitive Decline Risk': {
        'rs429358': {
            'gene': 'APOE',
            'variant': 'APOE ε4',
            'risk_allele': 'C',
            'description': 'APOE ε4 allele increases risk of cognitive decline and Alzheimer\'s',
            'interpretation': {
                'CC': 'ε4/ε4 - Highest risk of cognitive decline',
                'CT': 'ε4/ε3 - Increased risk',
                'TT': 'ε3/ε3 - Lower risk (most common)'
            },
            'cognitive_impact': 'Cognitive decline risk, Alzheimer\'s disease risk'
        },
        'rs7412': {
            'gene': 'APOE',
            'variant': 'APOE ε2',
            'risk_allele': 'T',
            'description': 'APOE ε2 allele may be protective against cognitive decline',
            'interpretation': {
                'TT': 'ε2/ε2 - Protective against cognitive decline',
                'CT': 'ε2/ε3 - Some protection',
                'CC': 'ε3/ε3 - Standard risk'
            },
            'cognitive_impact': 'Cognitive decline protection'
        },
    },

    # Processing Speed (if available)
    'Processing Speed': {
        'rs4851266': {
            'gene': 'Unknown',
            'variant': 'rs4851266-T',
            'risk_allele': 'T',
            'description': 'Associated with higher Verbal IQ and processing speed',
            'interpretation': {
                'TT': 'Higher verbal IQ and processing speed',
                'CT': 'Intermediate',
                'CC': 'Standard processing speed'
            },
            'cognitive_impact': 'Verbal IQ and processing speed'
        },
    },
}

# Cognitive domain weights for scoring
COGNITIVE_DOMAIN_WEIGHTS = {
    'Working Memory': 0.20,
    'Long-term Memory': 0.25,
    'Learning Ability': 0.25,
    'General Intelligence': 0.20,
    'Cognitive Decline Risk': 0.10,
}

# Learning style interpretations
LEARNING_STYLE_INTERPRETATIONS = {
    'Visual': 'May benefit from visual learning materials, diagrams, charts',
    'Auditory': 'May benefit from lectures, discussions, audio materials',
    'Kinesthetic': 'May benefit from hands-on activities, movement-based learning',
    'Mixed': 'Flexible learning style, adapts to different methods'
}
