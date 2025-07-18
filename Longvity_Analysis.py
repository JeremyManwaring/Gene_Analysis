#!/usr/bin/env python3
"""Longevity Analysis

This module merges the features from ``Longevity_Aging.py``,
``Longevity_Biomarkers.py`` and ``Probability_Calculation_Explanation.py``.
It exposes both the aging and biomarker analysers alongside helper
functions that demonstrate how the probability scores are derived.
"""

import sys
import subprocess
import os

# Auto-install pandas if missing
try:
    import pandas as pd
except ImportError:
    print("pandas not found. Installing pandas...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"])
    import pandas as pd

# Get genome file from command-line or default to genome.txt
if len(sys.argv) > 1:
    genome_file = sys.argv[1]
else:
    genome_file = "genome.txt"

if not os.path.isfile(genome_file):
    print(f"ERROR: Genome file '{genome_file}' not found. Please provide the correct file as an argument.")
    sys.exit(1)

import random
from typing import Dict, List


class LongevityAgingAnalyzer:
    """Analyze aging related genetic markers."""

    def __init__(self) -> None:
        # Telomere length and cellular aging markers
        self.telomere_markers = {
            'rs2736100': {'gene': 'TERT', 'effect': 'Telomere length', 'weight': 0.15},
            'rs7726159': {'gene': 'TERT', 'effect': 'Telomere length', 'weight': 0.12},
            'rs1317082': {'gene': 'TERC', 'effect': 'Telomere maintenance', 'weight': 0.10},
            'rs12696304': {'gene': 'TERC', 'effect': 'Telomere length', 'weight': 0.08},
            'rs10936599': {'gene': 'TERC', 'effect': 'Telomere maintenance', 'weight': 0.09},
            'rs16847897': {'gene': 'TERT', 'effect': 'Telomere length', 'weight': 0.11},
            'rs755017': {'gene': 'TERT', 'effect': 'Telomere maintenance', 'weight': 0.07},
            'rs4387287': {'gene': 'OBFC1', 'effect': 'Telomere length', 'weight': 0.06},
            'rs9420907': {'gene': 'OBFC1', 'effect': 'Telomere maintenance', 'weight': 0.05},
            'rs8105767': {'gene': 'ZNF208', 'effect': 'Telomere length', 'weight': 0.04},
        }

        # Oxidative stress markers
        self.oxidative_markers = {
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Oxidative stress', 'weight': 0.12},
            'rs1801131': {'gene': 'MTHFR', 'effect': 'Cellular damage', 'weight': 0.10},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Oxidative stress', 'weight': 0.08},
            'rs2070744': {'gene': 'NOS3', 'effect': 'Cellular damage', 'weight': 0.07},
            'rs1042522': {'gene': 'TP53', 'effect': 'Cellular repair', 'weight': 0.15},
            'rs1800371': {'gene': 'TP53', 'effect': 'DNA repair', 'weight': 0.13},
            'rs1045485': {'gene': 'CASP8', 'effect': 'Apoptosis', 'weight': 0.09},
            'rs3834129': {'gene': 'CASP8', 'effect': 'Cell death', 'weight': 0.08},
            'rs1800629': {'gene': 'TNF', 'effect': 'Inflammation', 'weight': 0.11},
            'rs361525': {'gene': 'TNF', 'effect': 'Inflammatory response', 'weight': 0.10},
        }

        # Inflammation markers
        self.inflammation_markers = {
            'rs1800795': {'gene': 'IL6', 'effect': 'Inflammation', 'weight': 0.14},
            'rs1800796': {'gene': 'IL6', 'effect': 'Immune aging', 'weight': 0.12},
            'rs1800896': {'gene': 'IL10', 'effect': 'Anti-inflammatory', 'weight': 0.10},
            'rs1800871': {'gene': 'IL10', 'effect': 'Immune regulation', 'weight': 0.09},
            'rs1800872': {'gene': 'IL10', 'effect': 'Inflammation control', 'weight': 0.08},
            'rs1800629': {'gene': 'TNF', 'effect': 'Pro-inflammatory', 'weight': 0.13},
            'rs361525': {'gene': 'TNF', 'effect': 'Inflammatory response', 'weight': 0.11},
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Inflammation', 'weight': 0.07},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Inflammatory response', 'weight': 0.06},
            'rs2070744': {'gene': 'NOS3', 'effect': 'Inflammation', 'weight': 0.05},
        }

        # Longevity markers
        self.longevity_markers = {
            'rs1042522': {'gene': 'TP53', 'effect': 'Lifespan regulation', 'weight': 0.18},
            'rs1800371': {'gene': 'TP53', 'effect': 'Longevity', 'weight': 0.16},
            'rs2736100': {'gene': 'TERT', 'effect': 'Cellular longevity', 'weight': 0.15},
            'rs7726159': {'gene': 'TERT', 'effect': 'Telomere longevity', 'weight': 0.13},
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Aging rate', 'weight': 0.12},
            'rs1801131': {'gene': 'MTHFR', 'effect': 'Lifespan', 'weight': 0.11},
            'rs1800795': {'gene': 'IL6', 'effect': 'Aging rate', 'weight': 0.10},
            'rs1800896': {'gene': 'IL10', 'effect': 'Longevity', 'weight': 0.09},
            'rs1045485': {'gene': 'CASP8', 'effect': 'Aging rate', 'weight': 0.08},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Lifespan', 'weight': 0.07},
        }

        # Cellular senescence markers
        self.senescence_markers = {
            'rs1042522': {'gene': 'TP53', 'effect': 'Senescence regulation', 'weight': 0.20},
            'rs1800371': {'gene': 'TP53', 'effect': 'Cellular senescence', 'weight': 0.18},
            'rs2736100': {'gene': 'TERT', 'effect': 'Senescence prevention', 'weight': 0.16},
            'rs7726159': {'gene': 'TERT', 'effect': 'Cellular aging', 'weight': 0.14},
            'rs1800795': {'gene': 'IL6', 'effect': 'Senescence', 'weight': 0.12},
            'rs1800896': {'gene': 'IL10', 'effect': 'Senescence control', 'weight': 0.10},
            'rs1045485': {'gene': 'CASP8', 'effect': 'Senescence', 'weight': 0.08},
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Cellular aging', 'weight': 0.07},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Senescence', 'weight': 0.06},
            'rs1800629': {'gene': 'TNF', 'effect': 'Cellular aging', 'weight': 0.05},
        }

        # Metabolic aging markers
        self.metabolic_markers = {
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Metabolic aging', 'weight': 0.15},
            'rs1801131': {'gene': 'MTHFR', 'effect': 'Metabolic rate', 'weight': 0.13},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Metabolic aging', 'weight': 0.12},
            'rs2070744': {'gene': 'NOS3', 'effect': 'Metabolic rate', 'weight': 0.11},
            'rs1800795': {'gene': 'IL6', 'effect': 'Metabolic aging', 'weight': 0.10},
            'rs1800896': {'gene': 'IL10', 'effect': 'Metabolic regulation', 'weight': 0.09},
            'rs1042522': {'gene': 'TP53', 'effect': 'Metabolic aging', 'weight': 0.08},
            'rs1800371': {'gene': 'TP53', 'effect': 'Metabolic rate', 'weight': 0.07},
            'rs1045485': {'gene': 'CASP8', 'effect': 'Metabolic aging', 'weight': 0.06},
            'rs1800629': {'gene': 'TNF', 'effect': 'Metabolic rate', 'weight': 0.05},
        }

    def load_genome_data(self, filename: str) -> Dict[str, str]:
        """Load genome data from a text file or generate random data."""
        genome_data: Dict[str, str] = {}
        try:
            with open(filename, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            genome_data[parts[0]] = parts[1]
        except FileNotFoundError:
            print(f"Warning: {filename} not found. Using simulated data.")
            all_snps = set()
            for markers in [
                self.telomere_markers,
                self.oxidative_markers,
                self.inflammation_markers,
                self.longevity_markers,
                self.senescence_markers,
                self.metabolic_markers,
            ]:
                all_snps.update(markers.keys())
            for snp in all_snps:
                genome_data[snp] = random.choice(['AA', 'AG', 'GG', 'AC', 'CC', 'AT', 'TT'])
        return genome_data

    @staticmethod
    def _effect_from_genotype(genotype: str) -> float:
        if 'A' in genotype and 'G' in genotype:
            return 0.5
        if genotype == 'GG':
            return 1.0
        if genotype == 'AA':
            return 0.0
        return 0.3

    def _analyze(self, genome_data: Dict[str, str], markers: Dict[str, Dict[str, float]], key: str) -> Dict:
        results = {
            'markers_found': 0,
            f'{key}_score': 0.0,
            'total_markers': len(markers),
            'risk_level': 'Unknown',
            'details': [],
        }
        for snp, info in markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                effect = self._effect_from_genotype(genome_data[snp])
                score = effect * info['weight']
                results[f'{key}_score'] += score
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genome_data[snp],
                    'effect': info['effect'],
                    'score': score,
                })
        score = results[f'{key}_score']
        if score >= 0.7:
            results['risk_level'] = 'Low'
        elif score >= 0.4:
            results['risk_level'] = 'Moderate'
        else:
            results['risk_level'] = 'High'
        return results

    def analyze_telomere_length(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.telomere_markers, 'telomere')

    def analyze_oxidative_stress(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.oxidative_markers, 'oxidative')

    def analyze_inflammation(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.inflammation_markers, 'inflammation')

    def analyze_longevity(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.longevity_markers, 'longevity')

    def analyze_senescence(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.senescence_markers, 'senescence')

    def analyze_metabolic_aging(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.metabolic_markers, 'metabolic')

    def calculate_polygenic_risk_score(self, genome_data: Dict[str, str]) -> Dict:
        """Calculate a simple polygenic risk score and probabilities."""
        all_markers = {}
        for collection in [
            self.telomere_markers,
            self.oxidative_markers,
            self.inflammation_markers,
            self.longevity_markers,
            self.senescence_markers,
            self.metabolic_markers,
        ]:
            all_markers.update(collection)

        total_score = 0.0
        markers_used = 0
        details = []
        for snp, info in all_markers.items():
            if snp in genome_data:
                markers_used += 1
                effect = self._effect_from_genotype(genome_data[snp])
                score = effect * info['weight']
                total_score += score
                details.append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genome_data[snp],
                    'effect': info['effect'],
                    'weight': info['weight'],
                    'score': score,
                })
        normalized_score = total_score / len(all_markers) if all_markers else 0
        aging_rate_prob = 1 - normalized_score
        longevity_prob = normalized_score
        cellular_health_prob = normalized_score
        return {
            'total_markers': len(all_markers),
            'markers_used': markers_used,
            'polygenic_score': normalized_score,
            'aging_rate_probability': aging_rate_prob,
            'longevity_probability': longevity_prob,
            'cellular_health_probability': cellular_health_prob,
            'overall_aging_risk': (
                'Low' if normalized_score >= 0.7 else
                'Moderate' if normalized_score >= 0.4 else 'High'
            ),
            'details': details,
        }

    def generate_recommendations(self, results: Dict) -> List[str]:
        """Create very basic lifestyle recommendations."""
        recs: List[str] = []
        if results['telomere']['risk_level'] == 'High':
            recs.extend([
                "Consider telomere-supporting supplements (vitamin D, omega-3)",
                "Implement stress reduction techniques (meditation, yoga)",
                "Regular exercise to support telomere maintenance",
                "Adequate sleep (7-9 hours per night)",
                "Anti-inflammatory diet rich in antioxidants",
            ])
        if results['oxidative']['risk_level'] == 'High':
            recs.extend([
                "High antioxidant diet (berries, dark chocolate, nuts)",
                "Consider antioxidant supplements (vitamin C, E, selenium)",
                "Reduce exposure to environmental toxins",
                "Regular detoxification practices",
                "Anti-inflammatory diet",
            ])
        if results['inflammation']['risk_level'] == 'High':
            recs.extend([
                "Anti-inflammatory diet (Mediterranean style)",
                "Omega-3 fatty acid supplementation",
                "Regular exercise to reduce inflammation",
                "Stress reduction techniques",
                "Adequate sleep and recovery",
            ])
        if results['longevity']['risk_level'] == 'High':
            recs.extend([
                "Caloric restriction or intermittent fasting",
                "Regular exercise (both cardio and strength training)",
                "Social connections and community engagement",
                "Continuous learning and mental stimulation",
                "Regular health monitoring and preventive care",
            ])
        if results['metabolic']['risk_level'] == 'High':
            recs.extend([
                "Regular exercise to maintain metabolic health",
                "Balanced diet with proper macronutrient ratios",
                "Adequate protein intake for muscle maintenance",
                "Regular strength training",
                "Monitor blood glucose and metabolic markers",
            ])
        if results['prs']['overall_aging_risk'] == 'High':
            recs.extend([
                "Comprehensive lifestyle modification program",
                "Regular medical check-ups and monitoring",
                "Consider working with a longevity specialist",
                "Implement all preventive measures",
                "Regular biomarker testing",
            ])
        return list(set(recs))

    def run_analysis(self, genome_file: str = "Genome.txt") -> Dict:
        data = self.load_genome_data(genome_file)
        results = {
            'telomere': self.analyze_telomere_length(data),
            'oxidative': self.analyze_oxidative_stress(data),
            'inflammation': self.analyze_inflammation(data),
            'longevity': self.analyze_longevity(data),
            'senescence': self.analyze_senescence(data),
            'metabolic': self.analyze_metabolic_aging(data),
            'prs': self.calculate_polygenic_risk_score(data),
        }
        return results


class LongevityBiomarkersAnalyzer:
    """Focus on positive longevity-promoting markers."""

    def __init__(self) -> None:
        self.longevity_markers = {
            'rs429358': {'gene': 'APOE', 'effect': 'APOE ε2 allele - longevity', 'weight': 0.20, 'protective_allele': 'T'},
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Folate metabolism - longevity', 'weight': 0.15, 'protective_allele': 'C'},
            'rs1801131': {'gene': 'MTHFR', 'effect': 'Homocysteine regulation', 'weight': 0.12, 'protective_allele': 'A'},
            'rs2736100': {'gene': 'TERT', 'effect': 'Telomere maintenance', 'weight': 0.18, 'protective_allele': 'G'},
            'rs7726159': {'gene': 'TERT', 'effect': 'Telomere length preservation', 'weight': 0.16, 'protective_allele': 'G'},
            'rs1317082': {'gene': 'TERC', 'effect': 'Telomere stability', 'weight': 0.14, 'protective_allele': 'G'},
            'rs1042522': {'gene': 'TP53', 'effect': 'Cellular repair - longevity', 'weight': 0.17, 'protective_allele': 'G'},
            'rs1800371': {'gene': 'TP53', 'effect': 'DNA repair enhancement', 'weight': 0.15, 'protective_allele': 'A'},
            'rs1800795': {'gene': 'IL6', 'effect': 'Anti-inflammatory - longevity', 'weight': 0.13, 'protective_allele': 'G'},
            'rs1800896': {'gene': 'IL10', 'effect': 'Immune regulation', 'weight': 0.11, 'protective_allele': 'A'},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Nitric oxide - longevity', 'weight': 0.12, 'protective_allele': 'G'},
            'rs2070744': {'gene': 'NOS3', 'effect': 'Vascular health', 'weight': 0.10, 'protective_allele': 'T'},
            'rs6265': {'gene': 'BDNF', 'effect': 'Brain health - longevity', 'weight': 0.14, 'protective_allele': 'G'},
            'rs1045485': {'gene': 'CASP8', 'effect': 'Apoptosis regulation', 'weight': 0.09, 'protective_allele': 'G'},
            'rs1800629': {'gene': 'TNF', 'effect': 'Inflammation control', 'weight': 0.08, 'protective_allele': 'A'},
        }

        self.centenarian_markers = {
            'rs429358': {'gene': 'APOE', 'effect': 'APOE ε2 - centenarian trait', 'weight': 0.25, 'protective_allele': 'T'},
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Methylation - longevity', 'weight': 0.20, 'protective_allele': 'C'},
            'rs2736100': {'gene': 'TERT', 'effect': 'Telomere - centenarian', 'weight': 0.22, 'protective_allele': 'G'},
            'rs1042522': {'gene': 'TP53', 'effect': 'Cellular repair - centenarian', 'weight': 0.21, 'protective_allele': 'G'},
            'rs1800795': {'gene': 'IL6', 'effect': 'Inflammation - centenarian', 'weight': 0.18, 'protective_allele': 'G'},
            'rs6265': {'gene': 'BDNF', 'effect': 'Cognitive longevity', 'weight': 0.19, 'protective_allele': 'G'},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Vascular longevity', 'weight': 0.17, 'protective_allele': 'G'},
            'rs1045485': {'gene': 'CASP8', 'effect': 'Cell survival - centenarian', 'weight': 0.15, 'protective_allele': 'G'},
            'rs1800629': {'gene': 'TNF', 'effect': 'Immune longevity', 'weight': 0.14, 'protective_allele': 'A'},
            'rs1800896': {'gene': 'IL10', 'effect': 'Anti-aging immune', 'weight': 0.16, 'protective_allele': 'A'},
        }

        self.anti_aging_markers = {
            'rs2736100': {'gene': 'TERT', 'effect': 'Telomere anti-aging', 'weight': 0.20, 'protective_allele': 'G'},
            'rs7726159': {'gene': 'TERT', 'effect': 'Cellular rejuvenation', 'weight': 0.18, 'protective_allele': 'G'},
            'rs1042522': {'gene': 'TP53', 'effect': 'Anti-cancer longevity', 'weight': 0.19, 'protective_allele': 'G'},
            'rs1800371': {'gene': 'TP53', 'effect': 'DNA repair anti-aging', 'weight': 0.17, 'protective_allele': 'A'},
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Methylation anti-aging', 'weight': 0.16, 'protective_allele': 'C'},
            'rs1800795': {'gene': 'IL6', 'effect': 'Anti-inflammatory aging', 'weight': 0.15, 'protective_allele': 'G'},
            'rs1800896': {'gene': 'IL10', 'effect': 'Immune anti-aging', 'weight': 0.14, 'protective_allele': 'A'},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Vascular anti-aging', 'weight': 0.13, 'protective_allele': 'G'},
            'rs6265': {'gene': 'BDNF', 'effect': 'Brain anti-aging', 'weight': 0.15, 'protective_allele': 'G'},
            'rs1045485': {'gene': 'CASP8', 'effect': 'Cell survival anti-aging', 'weight': 0.12, 'protective_allele': 'G'},
        }

        self.metabolic_longevity_markers = {
            'rs1801133': {'gene': 'MTHFR', 'effect': 'Metabolic longevity', 'weight': 0.18, 'protective_allele': 'C'},
            'rs1799983': {'gene': 'NOS3', 'effect': 'Metabolic health', 'weight': 0.16, 'protective_allele': 'G'},
            'rs1800795': {'gene': 'IL6', 'effect': 'Metabolic anti-aging', 'weight': 0.15, 'protective_allele': 'G'},
            'rs1800896': {'gene': 'IL10', 'effect': 'Metabolic regulation', 'weight': 0.14, 'protective_allele': 'A'},
            'rs1042522': {'gene': 'TP53', 'effect': 'Metabolic repair', 'weight': 0.17, 'protective_allele': 'G'},
            'rs6265': {'gene': 'BDNF', 'effect': 'Metabolic brain health', 'weight': 0.13, 'protective_allele': 'G'},
            'rs2736100': {'gene': 'TERT', 'effect': 'Metabolic telomere', 'weight': 0.16, 'protective_allele': 'G'},
            'rs1045485': {'gene': 'CASP8', 'effect': 'Metabolic survival', 'weight': 0.11, 'protective_allele': 'G'},
            'rs1800629': {'gene': 'TNF', 'effect': 'Metabolic inflammation', 'weight': 0.12, 'protective_allele': 'A'},
            'rs1800371': {'gene': 'TP53', 'effect': 'Metabolic DNA repair', 'weight': 0.14, 'protective_allele': 'A'},
        }

    def load_genome_data(self, filename: str) -> Dict[str, str]:
        genome_data: Dict[str, str] = {}
        try:
            with open(filename, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            genome_data[parts[0]] = parts[1]
        except FileNotFoundError:
            print(f"Warning: {filename} not found. Using simulated data.")
            all_snps = set()
            for markers in [
                self.longevity_markers,
                self.centenarian_markers,
                self.anti_aging_markers,
                self.metabolic_longevity_markers,
            ]:
                all_snps.update(markers.keys())
            for snp in all_snps:
                genome_data[snp] = random.choice(['AA', 'AG', 'GG', 'AC', 'CC', 'AT', 'TT', 'CT', 'TC'])
        return genome_data

    @staticmethod
    def _effect(genotype: str, allele: str) -> float:
        if allele in genotype:
            return 1.0 if genotype.count(allele) == 2 else 0.7
        return 0.0

    def _analyze(self, genome_data: Dict[str, str], markers: Dict[str, Dict], key: str) -> Dict:
        results = {
            'markers_found': 0,
            'total_markers': len(markers),
            f'{key}_score': 0.0,
            'protective_alleles_found': 0,
            f'{key}_level': 'Unknown',
            'details': [],
        }
        for snp, info in markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                has_p = info['protective_allele'] in genotype
                if has_p:
                    results['protective_alleles_found'] += 1
                effect = self._effect(genotype, info['protective_allele'])
                score = effect * info['weight']
                results[f'{key}_score'] += score
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'protective_allele': info['protective_allele'],
                    'has_protective': has_p,
                    'score': score,
                })
        score = results[f'{key}_score']
        if score >= 0.8:
            level = 'Exceptional'
        elif score >= 0.6:
            level = 'High'
        elif score >= 0.4:
            level = 'Moderate'
        elif score >= 0.2:
            level = 'Low'
        else:
            level = 'Poor'
        results[f'{key}_level'] = level
        return results

    def analyze_longevity_biomarkers(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.longevity_markers, 'longevity')

    def analyze_centenarian_traits(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.centenarian_markers, 'centenarian')

    def analyze_anti_aging_markers(self, data: Dict[str, str]) -> Dict:
        return self._analyze(data, self.anti_aging_markers, 'anti_aging')

    def calculate_longevity_probability(self, genome_data: Dict[str, str]) -> Dict:
        all_markers = {}
        for markers in [
            self.longevity_markers,
            self.centenarian_markers,
            self.anti_aging_markers,
            self.metabolic_longevity_markers,
        ]:
            all_markers.update(markers)
        total_score = 0.0
        markers_used = 0
        protective = 0
        details = []
        for snp, info in all_markers.items():
            if snp in genome_data:
                markers_used += 1
                genotype = genome_data[snp]
                has_p = info['protective_allele'] in genotype
                if has_p:
                    protective += 1
                effect = self._effect(genotype, info['protective_allele'])
                score = effect * info['weight']
                total_score += score
                details.append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'protective_allele': info['protective_allele'],
                    'has_protective': has_p,
                    'score': score,
                })
        normalized_score = total_score / len(all_markers) if all_markers else 0
        longevity_probability = normalized_score
        centenarian_probability = normalized_score * 0.8
        anti_aging_probability = normalized_score
        healthy_aging_probability = normalized_score * 0.9
        return {
            'total_markers': len(all_markers),
            'markers_used': markers_used,
            'protective_alleles_found': protective,
            'longevity_score': normalized_score,
            'longevity_probability': longevity_probability,
            'centenarian_probability': centenarian_probability,
            'anti_aging_probability': anti_aging_probability,
            'healthy_aging_probability': healthy_aging_probability,
            'overall_longevity_potential': (
                'Exceptional' if normalized_score >= 0.8 else
                'High' if normalized_score >= 0.6 else
                'Moderate' if normalized_score >= 0.4 else
                'Low' if normalized_score >= 0.2 else 'Poor'
            ),
            'details': details,
        }

    def run_analysis(self, genome_file: str = "Genome.txt") -> Dict:
        data = self.load_genome_data(genome_file)
        results = {
            'longevity': self.analyze_longevity_biomarkers(data),
            'centenarian': self.analyze_centenarian_traits(data),
            'anti_aging': self.analyze_anti_aging_markers(data),
            'probability': self.calculate_longevity_probability(data),
        }
        return results


def explain_probability_calculation() -> Dict:
    """Replicate the probability calculation explanation from the old script."""
    print("🔬 PROBABILITY CALCULATION BREAKDOWN")
    print("=" * 60)
    print("\n📊 STEP 1: GENOTYPE EFFECT CALCULATION")
    print("-" * 40)
    print("For each SNP, effects are calculated based on genotype:")
    print("• GG (Protective): effect = 1.0")
    print("• AG/GA (Heterozygous): effect = 0.5")
    print("• AA (Risk): effect = 0.0")
    print("• Other genotypes: effect = 0.3")

    example_genotypes = ['GG', 'AG', 'AA', 'CT']
    example_weights = [0.15, 0.12, 0.10, 0.08]
    print("\n📋 EXAMPLE CALCULATIONS:")
    for i, (genotype, weight) in enumerate(zip(example_genotypes, example_weights)):
        if 'A' in genotype and 'G' in genotype:
            effect = 0.5
        elif genotype == 'GG':
            effect = 1.0
        elif genotype == 'AA':
            effect = 0.0
        else:
            effect = 0.3
        score = effect * weight
        print(f"SNP {i+1}: {genotype} × {weight} = {score:.3f}")

    print("\n📊 STEP 2: POLYGENIC SCORE CALCULATION")
    print("-" * 40)
    total_markers = 25
    markers_used = 21
    polygenic_score = 0.024
    print(f"Total markers available: {total_markers}")
    print(f"Markers found in genome: {markers_used}")
    print("Individual SNP scores summed up")
    print(f"Normalized by total markers: {polygenic_score:.3f}")

    print("\n📊 STEP 3: PROBABILITY CALCULATIONS")
    print("-" * 40)
    aging_rate_prob = 1 - polygenic_score
    longevity_prob = polygenic_score
    cellular_health_prob = polygenic_score
    print(f"Formula 1: Aging Rate Probability = 1 - {polygenic_score:.3f}")
    print(f"         = {aging_rate_prob:.1%}")
    print(f"\nFormula 2: Longevity Probability = {longevity_prob:.3f}")
    print(f"         = {longevity_prob:.1%}")
    print(f"\nFormula 3: Cellular Health Probability = {cellular_health_prob:.3f}")
    print(f"         = {cellular_health_prob:.1%}")

    print("\n📊 STEP 4: INTERPRETATION")
    print("-" * 40)
    print("• Higher polygenic score = better aging genetics")
    print(f"• Your score: {polygenic_score:.3f} (very low)")
    print(f"• This means:")
    print(f"  - High risk of accelerated aging ({aging_rate_prob:.1%})")
    print(f"  - Low probability of extended lifespan ({longevity_prob:.1%})")
    print(f"  - Low cellular health maintenance ({cellular_health_prob:.1%})")
    return {
        'polygenic_score': polygenic_score,
        'aging_rate_probability': aging_rate_prob,
        'longevity_probability': longevity_prob,
        'cellular_health_probability': cellular_health_prob,
    }


def show_detailed_snp_breakdown() -> None:
    """Demonstrate SNP by SNP scoring."""
    print("\n🔬 DETAILED SNP BREAKDOWN")
    print("=" * 60)
    example_snps = [
        {'snp': 'rs2736100', 'gene': 'TERT', 'genotype': 'AA', 'weight': 0.15},
        {'snp': 'rs7726159', 'gene': 'TERT', 'genotype': 'AG', 'weight': 0.12},
        {'snp': 'rs1317082', 'gene': 'TERC', 'genotype': 'GG', 'weight': 0.10},
        {'snp': 'rs12696304', 'gene': 'TERC', 'genotype': 'AA', 'weight': 0.08},
        {'snp': 'rs10936599', 'gene': 'TERC', 'genotype': 'CT', 'weight': 0.09},
    ]
    total_score = 0.0
    for snp in example_snps:
        genotype = snp['genotype']
        if 'A' in genotype and 'G' in genotype:
            effect = 0.5
        elif genotype == 'GG':
            effect = 1.0
        elif genotype == 'AA':
            effect = 0.0
        else:
            effect = 0.3
        score = effect * snp['weight']
        total_score += score
        print(f"{snp['snp']} ({snp['gene']}): {genotype} -> effect={effect} -> score={score:.3f}")
    print(f"\nSubtotal for these 5 SNPs: {total_score:.3f}")
    print("With 25 total SNPs, normalized score would be very low")
    print("This explains the 97.6% aging rate probability")


def main() -> None:
    aging_analyzer = LongevityAgingAnalyzer()
    aging_results = aging_analyzer.run_analysis()
    print("\nAging analysis complete. Overall risk:", aging_results['prs']['overall_aging_risk'])

    biomarkers_analyzer = LongevityBiomarkersAnalyzer()
    biomarker_results = biomarkers_analyzer.run_analysis()
    print("\nBiomarker analysis complete. Overall potential:", biomarker_results['probability']['overall_longevity_potential'])

    explain_probability_calculation()
    show_detailed_snp_breakdown()


if __name__ == "__main__":
    main()
