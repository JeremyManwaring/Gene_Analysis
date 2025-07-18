#!/usr/bin/env python3
"""
Comprehensive Longevity and Aging Analysis
Analyzes genetic markers related to aging, longevity, and cellular health
"""

import random
import math
from typing import Dict, List, Tuple

class LongevityAgingAnalyzer:
    def __init__(self):
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
            'rs8105767': {'gene': 'ZNF208', 'effect': 'Telomere length', 'weight': 0.04}
        }
        
        # Oxidative stress and cellular damage markers
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
            'rs361525': {'gene': 'TNF', 'effect': 'Inflammatory response', 'weight': 0.10}
        }
        
        # Inflammation and immune aging markers
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
            'rs2070744': {'gene': 'NOS3', 'effect': 'Inflammation', 'weight': 0.05}
        }
        
        # Longevity and lifespan markers
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
            'rs2802292': {'gene': 'FOXO3', 'effect': 'Longevity variant', 'weight': 0.15}
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
            'rs2802292': {'gene': 'FOXO3', 'effect': 'Longevity variant', 'weight': 0.14}
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
            'rs1800629': {'gene': 'TNF', 'effect': 'Metabolic rate', 'weight': 0.05}
            ,'rs2802292': {'gene': 'FOXO3', 'effect': 'Metabolic longevity', 'weight': 0.13}
        }

    def load_genome_data(self, filename: str) -> Dict[str, str]:
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
            all_snps = set()
            for markers in [self.telomere_markers, self.oxidative_markers, 
                          self.inflammation_markers, self.longevity_markers,
                          self.senescence_markers, self.metabolic_markers]:
                all_snps.update(markers.keys())
            
            for snp in all_snps:
                genotypes = ['AA', 'AG', 'GG', 'AC', 'CC', 'AT', 'TT']
                genome_data[snp] = random.choice(genotypes)
        
        return genome_data

    def analyze_telomere_length(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze telomere length markers"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.telomere_markers),
            'telomere_score': 0.0,
            'risk_level': 'Unknown',
            'details': []
        }
        
        for snp, info in self.telomere_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Calculate effect based on genotype
                if 'A' in genotype and 'G' in genotype:  # Heterozygous
                    effect = 0.5
                elif genotype == 'GG':  # Protective
                    effect = 1.0
                elif genotype == 'AA':  # Risk
                    effect = 0.0
                else:
                    effect = 0.3  # Other genotypes
                
                score = effect * info['weight']
                results['telomere_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'score': score
                })
        
        # Determine risk level
        if results['telomere_score'] >= 0.7:
            results['risk_level'] = 'Low'
        elif results['telomere_score'] >= 0.4:
            results['risk_level'] = 'Moderate'
        else:
            results['risk_level'] = 'High'
        
        return results

    def analyze_oxidative_stress(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze oxidative stress markers"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.oxidative_markers),
            'oxidative_score': 0.0,
            'risk_level': 'Unknown',
            'details': []
        }
        
        for snp, info in self.oxidative_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Calculate effect based on genotype
                if 'A' in genotype and 'G' in genotype:  # Heterozygous
                    effect = 0.5
                elif genotype == 'GG':  # Protective
                    effect = 1.0
                elif genotype == 'AA':  # Risk
                    effect = 0.0
                else:
                    effect = 0.3  # Other genotypes
                
                score = effect * info['weight']
                results['oxidative_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'score': score
                })
        
        # Determine risk level
        if results['oxidative_score'] >= 0.7:
            results['risk_level'] = 'Low'
        elif results['oxidative_score'] >= 0.4:
            results['risk_level'] = 'Moderate'
        else:
            results['risk_level'] = 'High'
        
        return results

    def analyze_inflammation(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze inflammation markers"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.inflammation_markers),
            'inflammation_score': 0.0,
            'risk_level': 'Unknown',
            'details': []
        }
        
        for snp, info in self.inflammation_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Calculate effect based on genotype
                if 'A' in genotype and 'G' in genotype:  # Heterozygous
                    effect = 0.5
                elif genotype == 'GG':  # Protective
                    effect = 1.0
                elif genotype == 'AA':  # Risk
                    effect = 0.0
                else:
                    effect = 0.3  # Other genotypes
                
                score = effect * info['weight']
                results['inflammation_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'score': score
                })
        
        # Determine risk level
        if results['inflammation_score'] >= 0.7:
            results['risk_level'] = 'Low'
        elif results['inflammation_score'] >= 0.4:
            results['risk_level'] = 'Moderate'
        else:
            results['risk_level'] = 'High'
        
        return results

    def analyze_longevity(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze longevity markers"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.longevity_markers),
            'longevity_score': 0.0,
            'risk_level': 'Unknown',
            'details': []
        }
        
        for snp, info in self.longevity_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Calculate effect based on genotype
                if 'A' in genotype and 'G' in genotype:  # Heterozygous
                    effect = 0.5
                elif genotype == 'GG':  # Protective
                    effect = 1.0
                elif genotype == 'AA':  # Risk
                    effect = 0.0
                else:
                    effect = 0.3  # Other genotypes
                
                score = effect * info['weight']
                results['longevity_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'score': score
                })
        
        # Determine risk level
        if results['longevity_score'] >= 0.7:
            results['risk_level'] = 'Low'
        elif results['longevity_score'] >= 0.4:
            results['risk_level'] = 'Moderate'
        else:
            results['risk_level'] = 'High'
        
        return results

    def analyze_senescence(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze cellular senescence markers"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.senescence_markers),
            'senescence_score': 0.0,
            'risk_level': 'Unknown',
            'details': []
        }
        
        for snp, info in self.senescence_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Calculate effect based on genotype
                if 'A' in genotype and 'G' in genotype:  # Heterozygous
                    effect = 0.5
                elif genotype == 'GG':  # Protective
                    effect = 1.0
                elif genotype == 'AA':  # Risk
                    effect = 0.0
                else:
                    effect = 0.3  # Other genotypes
                
                score = effect * info['weight']
                results['senescence_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'score': score
                })
        
        # Determine risk level
        if results['senescence_score'] >= 0.7:
            results['risk_level'] = 'Low'
        elif results['senescence_score'] >= 0.4:
            results['risk_level'] = 'Moderate'
        else:
            results['risk_level'] = 'High'
        
        return results

    def analyze_metabolic_aging(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze metabolic aging markers"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.metabolic_markers),
            'metabolic_score': 0.0,
            'risk_level': 'Unknown',
            'details': []
        }
        
        for snp, info in self.metabolic_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Calculate effect based on genotype
                if 'A' in genotype and 'G' in genotype:  # Heterozygous
                    effect = 0.5
                elif genotype == 'GG':  # Protective
                    effect = 1.0
                elif genotype == 'AA':  # Risk
                    effect = 0.0
                else:
                    effect = 0.3  # Other genotypes
                
                score = effect * info['weight']
                results['metabolic_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'score': score
                })
        
        # Determine risk level
        if results['metabolic_score'] >= 0.7:
            results['risk_level'] = 'Low'
        elif results['metabolic_score'] >= 0.4:
            results['risk_level'] = 'Moderate'
        else:
            results['risk_level'] = 'High'
        
        return results

    def calculate_polygenic_risk_score(self, genome_data: Dict[str, str]) -> Dict:
        """Calculate comprehensive polygenic risk score for aging"""
        all_markers = {}
        all_markers.update(self.telomere_markers)
        all_markers.update(self.oxidative_markers)
        all_markers.update(self.inflammation_markers)
        all_markers.update(self.longevity_markers)
        all_markers.update(self.senescence_markers)
        all_markers.update(self.metabolic_markers)
        
        total_score = 0.0
        markers_used = 0
        details = []
        
        for snp, info in all_markers.items():
            if snp in genome_data:
                markers_used += 1
                genotype = genome_data[snp]
                
                # Calculate effect based on genotype
                if 'A' in genotype and 'G' in genotype:  # Heterozygous
                    effect = 0.5
                elif genotype == 'GG':  # Protective
                    effect = 1.0
                elif genotype == 'AA':  # Risk
                    effect = 0.0
                else:
                    effect = 0.3  # Other genotypes
                
                score = effect * info['weight']
                total_score += score
                
                details.append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'weight': info['weight'],
                    'score': score
                })
        
        # Normalize score to 0-1 range
        normalized_score = total_score / len(all_markers) if all_markers else 0
        
        # Calculate probability scores
        aging_rate_prob = 1 - normalized_score  # Higher score = lower aging rate
        longevity_prob = normalized_score  # Higher score = higher longevity
        cellular_health_prob = normalized_score  # Higher score = better cellular health
        
        return {
            'total_markers': len(all_markers),
            'markers_used': markers_used,
            'polygenic_score': normalized_score,
            'aging_rate_probability': aging_rate_prob,
            'longevity_probability': longevity_prob,
            'cellular_health_probability': cellular_health_prob,
            'overall_aging_risk': 'Low' if normalized_score >= 0.7 else 'Moderate' if normalized_score >= 0.4 else 'High',
            'details': details
        }

    def generate_recommendations(self, results: Dict) -> List[str]:
        """Generate personalized recommendations based on analysis"""
        recommendations = []
        
        # Telomere recommendations
        if results['telomere']['risk_level'] == 'High':
            recommendations.extend([
                "Consider telomere-supporting supplements (vitamin D, omega-3)",
                "Implement stress reduction techniques (meditation, yoga)",
                "Regular exercise to support telomere maintenance",
                "Adequate sleep (7-9 hours per night)",
                "Anti-inflammatory diet rich in antioxidants"
            ])
        elif results['telomere']['risk_level'] == 'Moderate':
            recommendations.extend([
                "Maintain healthy lifestyle habits",
                "Regular moderate exercise",
                "Stress management techniques",
                "Balanced diet with antioxidants"
            ])
        
        # Oxidative stress recommendations
        if results['oxidative']['risk_level'] == 'High':
            recommendations.extend([
                "High antioxidant diet (berries, dark chocolate, nuts)",
                "Consider antioxidant supplements (vitamin C, E, selenium)",
                "Reduce exposure to environmental toxins",
                "Regular detoxification practices",
                "Anti-inflammatory diet"
            ])
        
        # Inflammation recommendations
        if results['inflammation']['risk_level'] == 'High':
            recommendations.extend([
                "Anti-inflammatory diet (Mediterranean style)",
                "Omega-3 fatty acid supplementation",
                "Regular exercise to reduce inflammation",
                "Stress reduction techniques",
                "Adequate sleep and recovery"
            ])
        
        # Longevity recommendations
        if results['longevity']['risk_level'] == 'High':
            recommendations.extend([
                "Caloric restriction or intermittent fasting",
                "Regular exercise (both cardio and strength training)",
                "Social connections and community engagement",
                "Continuous learning and mental stimulation",
                "Regular health monitoring and preventive care"
            ])
        
        # Metabolic aging recommendations
        if results['metabolic']['risk_level'] == 'High':
            recommendations.extend([
                "Regular exercise to maintain metabolic health",
                "Balanced diet with proper macronutrient ratios",
                "Adequate protein intake for muscle maintenance",
                "Regular strength training",
                "Monitor blood glucose and metabolic markers"
            ])
        
        # General recommendations based on overall score
        if results['prs']['overall_aging_risk'] == 'High':
            recommendations.extend([
                "Comprehensive lifestyle modification program",
                "Regular medical check-ups and monitoring",
                "Consider working with a longevity specialist",
                "Implement all preventive measures",
                "Regular biomarker testing"
            ])
        
        return list(set(recommendations))  # Remove duplicates

    def run_analysis(self, genome_file: str = "Genome.txt") -> Dict:
        """Run comprehensive longevity and aging analysis"""
        print("🧬 LONGEVITY & AGING ANALYSIS")
        print("=" * 50)
        
        # Load genome data
        genome_data = self.load_genome_data(genome_file)
        
        # Run all analyses
        telomere_results = self.analyze_telomere_length(genome_data)
        oxidative_results = self.analyze_oxidative_stress(genome_data)
        inflammation_results = self.analyze_inflammation(genome_data)
        longevity_results = self.analyze_longevity(genome_data)
        senescence_results = self.analyze_senescence(genome_data)
        metabolic_results = self.analyze_metabolic_aging(genome_data)
        prs_results = self.calculate_polygenic_risk_score(genome_data)
        
        # Compile results
        results = {
            'telomere': telomere_results,
            'oxidative': oxidative_results,
            'inflammation': inflammation_results,
            'longevity': longevity_results,
            'senescence': senescence_results,
            'metabolic': metabolic_results,
            'prs': prs_results
        }
        
        # Print detailed results
        print(f"\n📊 TELOMERE LENGTH ANALYSIS")
        print(f"Markers found: {telomere_results['markers_found']}/{telomere_results['total_markers']}")
        print(f"Telomere Score: {telomere_results['telomere_score']:.3f}")
        print(f"Risk Level: {telomere_results['risk_level']}")
        
        print(f"\n📊 OXIDATIVE STRESS ANALYSIS")
        print(f"Markers found: {oxidative_results['markers_found']}/{oxidative_results['total_markers']}")
        print(f"Oxidative Score: {oxidative_results['oxidative_score']:.3f}")
        print(f"Risk Level: {oxidative_results['risk_level']}")
        
        print(f"\n📊 INFLAMMATION ANALYSIS")
        print(f"Markers found: {inflammation_results['markers_found']}/{inflammation_results['total_markers']}")
        print(f"Inflammation Score: {inflammation_results['inflammation_score']:.3f}")
        print(f"Risk Level: {inflammation_results['risk_level']}")
        
        print(f"\n📊 LONGEVITY ANALYSIS")
        print(f"Markers found: {longevity_results['markers_found']}/{longevity_results['total_markers']}")
        print(f"Longevity Score: {longevity_results['longevity_score']:.3f}")
        print(f"Risk Level: {longevity_results['risk_level']}")
        
        print(f"\n📊 CELLULAR SENESCENCE ANALYSIS")
        print(f"Markers found: {senescence_results['markers_found']}/{senescence_results['total_markers']}")
        print(f"Senescence Score: {senescence_results['senescence_score']:.3f}")
        print(f"Risk Level: {senescence_results['risk_level']}")
        
        print(f"\n📊 METABOLIC AGING ANALYSIS")
        print(f"Markers found: {metabolic_results['markers_found']}/{metabolic_results['total_markers']}")
        print(f"Metabolic Score: {metabolic_results['metabolic_score']:.3f}")
        print(f"Risk Level: {metabolic_results['risk_level']}")
        
        print(f"\n🎯 POLYGENIC RISK SCORE SUMMARY")
        print(f"Total markers analyzed: {prs_results['total_markers']}")
        print(f"Markers used: {prs_results['markers_used']}")
        print(f"Polygenic Score: {prs_results['polygenic_score']:.3f}")
        print(f"Overall Aging Risk: {prs_results['overall_aging_risk']}")
        
        print(f"\n📈 PROBABILITY REPORT")
        print(f"Aging Rate Probability: {prs_results['aging_rate_probability']:.1%}")
        print(f"Longevity Probability: {prs_results['longevity_probability']:.1%}")
        print(f"Cellular Health Probability: {prs_results['cellular_health_probability']:.1%}")
        
        # Generate and print recommendations
        recommendations = self.generate_recommendations(results)
        
        print(f"\n💡 PERSONALIZED RECOMMENDATIONS")
        for i, rec in enumerate(recommendations, 1):
            print(f"{i}. {rec}")
        
        return results

def main():
    """Main function to run the analysis"""
    analyzer = LongevityAgingAnalyzer()
    results = analyzer.run_analysis()
    
    print(f"\n✅ Analysis complete!")
    print(f"📋 Summary: {results['prs']['overall_aging_risk']} overall aging risk")
    print(f"🎯 Key focus areas identified based on genetic markers")

if __name__ == "__main__":
    main() 