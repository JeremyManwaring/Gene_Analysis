#!/usr/bin/env python3
"""
Longevity Biomarkers Analysis
Focuses on positive longevity-promoting markers and protective factors
"""

import random
from typing import Dict, List

class LongevityBiomarkersAnalyzer:
    def __init__(self):
        # Longevity-promoting SNPs (positive markers)
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
            'rs1800629': {'gene': 'TNF', 'effect': 'Inflammation control', 'weight': 0.08, 'protective_allele': 'A'}
        }
        
        # Centenarian-associated SNPs
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
            'rs1800896': {'gene': 'IL10', 'effect': 'Anti-aging immune', 'weight': 0.16, 'protective_allele': 'A'}
        }
        
        # Anti-aging and cellular health markers
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
            'rs1045485': {'gene': 'CASP8', 'effect': 'Cell survival anti-aging', 'weight': 0.12, 'protective_allele': 'G'}
        }
        
        # Metabolic health and longevity markers
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
            'rs1800371': {'gene': 'TP53', 'effect': 'Metabolic DNA repair', 'weight': 0.14, 'protective_allele': 'A'}
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
            for markers in [self.longevity_markers, self.centenarian_markers, 
                          self.anti_aging_markers, self.metabolic_longevity_markers]:
                all_snps.update(markers.keys())
            
            for snp in all_snps:
                genotypes = ['AA', 'AG', 'GG', 'AC', 'CC', 'AT', 'TT', 'CT', 'TC']
                genome_data[snp] = random.choice(genotypes)
        
        return genome_data

    def analyze_longevity_biomarkers(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze positive longevity biomarkers"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.longevity_markers),
            'longevity_score': 0.0,
            'protective_alleles_found': 0,
            'longevity_level': 'Unknown',
            'details': []
        }
        
        for snp, info in self.longevity_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Check for protective alleles
                has_protective = info['protective_allele'] in genotype
                if has_protective:
                    results['protective_alleles_found'] += 1
                
                # Calculate effect based on genotype
                if info['protective_allele'] in genotype:
                    if genotype.count(info['protective_allele']) == 2:  # Homozygous protective
                        effect = 1.0
                    else:  # Heterozygous protective
                        effect = 0.7
                else:  # No protective alleles
                    effect = 0.0
                
                score = effect * info['weight']
                results['longevity_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'protective_allele': info['protective_allele'],
                    'has_protective': has_protective,
                    'score': score
                })
        
        # Determine longevity level
        if results['longevity_score'] >= 0.8:
            results['longevity_level'] = 'Exceptional'
        elif results['longevity_score'] >= 0.6:
            results['longevity_level'] = 'High'
        elif results['longevity_score'] >= 0.4:
            results['longevity_level'] = 'Moderate'
        elif results['longevity_score'] >= 0.2:
            results['longevity_level'] = 'Low'
        else:
            results['longevity_level'] = 'Poor'
        
        return results

    def analyze_centenarian_traits(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze centenarian-associated traits"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.centenarian_markers),
            'centenarian_score': 0.0,
            'protective_alleles_found': 0,
            'centenarian_potential': 'Unknown',
            'details': []
        }
        
        for snp, info in self.centenarian_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Check for protective alleles
                has_protective = info['protective_allele'] in genotype
                if has_protective:
                    results['protective_alleles_found'] += 1
                
                # Calculate effect based on genotype
                if info['protective_allele'] in genotype:
                    if genotype.count(info['protective_allele']) == 2:  # Homozygous protective
                        effect = 1.0
                    else:  # Heterozygous protective
                        effect = 0.7
                else:  # No protective alleles
                    effect = 0.0
                
                score = effect * info['weight']
                results['centenarian_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'protective_allele': info['protective_allele'],
                    'has_protective': has_protective,
                    'score': score
                })
        
        # Determine centenarian potential
        if results['centenarian_score'] >= 0.8:
            results['centenarian_potential'] = 'Exceptional'
        elif results['centenarian_score'] >= 0.6:
            results['centenarian_potential'] = 'High'
        elif results['centenarian_score'] >= 0.4:
            results['centenarian_potential'] = 'Moderate'
        elif results['centenarian_score'] >= 0.2:
            results['centenarian_potential'] = 'Low'
        else:
            results['centenarian_potential'] = 'Poor'
        
        return results

    def analyze_anti_aging_markers(self, genome_data: Dict[str, str]) -> Dict:
        """Analyze anti-aging markers"""
        results = {
            'markers_found': 0,
            'total_markers': len(self.anti_aging_markers),
            'anti_aging_score': 0.0,
            'protective_alleles_found': 0,
            'anti_aging_potential': 'Unknown',
            'details': []
        }
        
        for snp, info in self.anti_aging_markers.items():
            if snp in genome_data:
                results['markers_found'] += 1
                genotype = genome_data[snp]
                
                # Check for protective alleles
                has_protective = info['protective_allele'] in genotype
                if has_protective:
                    results['protective_alleles_found'] += 1
                
                # Calculate effect based on genotype
                if info['protective_allele'] in genotype:
                    if genotype.count(info['protective_allele']) == 2:  # Homozygous protective
                        effect = 1.0
                    else:  # Heterozygous protective
                        effect = 0.7
                else:  # No protective alleles
                    effect = 0.0
                
                score = effect * info['weight']
                results['anti_aging_score'] += score
                
                results['details'].append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'effect': info['effect'],
                    'protective_allele': info['protective_allele'],
                    'has_protective': has_protective,
                    'score': score
                })
        
        # Determine anti-aging potential
        if results['anti_aging_score'] >= 0.8:
            results['anti_aging_potential'] = 'Exceptional'
        elif results['anti_aging_score'] >= 0.6:
            results['anti_aging_potential'] = 'High'
        elif results['anti_aging_score'] >= 0.4:
            results['anti_aging_potential'] = 'Moderate'
        elif results['anti_aging_score'] >= 0.2:
            results['anti_aging_potential'] = 'Low'
        else:
            results['anti_aging_potential'] = 'Poor'
        
        return results

    def calculate_longevity_probability(self, genome_data: Dict[str, str]) -> Dict:
        """Calculate positive longevity probability scores"""
        all_markers = {}
        all_markers.update(self.longevity_markers)
        all_markers.update(self.centenarian_markers)
        all_markers.update(self.anti_aging_markers)
        all_markers.update(self.metabolic_longevity_markers)
        
        total_score = 0.0
        markers_used = 0
        protective_alleles = 0
        details = []
        
        for snp, info in all_markers.items():
            if snp in genome_data:
                markers_used += 1
                genotype = genome_data[snp]
                
                # Check for protective alleles
                has_protective = info['protective_allele'] in genotype
                if has_protective:
                    protective_alleles += 1
                
                # Calculate effect based on genotype
                if info['protective_allele'] in genotype:
                    if genotype.count(info['protective_allele']) == 2:  # Homozygous protective
                        effect = 1.0
                    else:  # Heterozygous protective
                        effect = 0.7
                else:  # No protective alleles
                    effect = 0.0
                
                score = effect * info['weight']
                total_score += score
                
                details.append({
                    'snp': snp,
                    'gene': info['gene'],
                    'genotype': genotype,
                    'protective_allele': info['protective_allele'],
                    'has_protective': has_protective,
                    'score': score
                })
        
        # Normalize score to 0-1 range
        normalized_score = total_score / len(all_markers) if all_markers else 0
        
        # Calculate positive probability scores
        longevity_probability = normalized_score  # Higher score = higher longevity
        centenarian_probability = normalized_score * 0.8  # Slightly lower for centenarian
        anti_aging_probability = normalized_score  # Higher score = better anti-aging
        healthy_aging_probability = normalized_score * 0.9  # Slightly lower for healthy aging
        
        return {
            'total_markers': len(all_markers),
            'markers_used': markers_used,
            'protective_alleles_found': protective_alleles,
            'longevity_score': normalized_score,
            'longevity_probability': longevity_probability,
            'centenarian_probability': centenarian_probability,
            'anti_aging_probability': anti_aging_probability,
            'healthy_aging_probability': healthy_aging_probability,
            'overall_longevity_potential': 'Exceptional' if normalized_score >= 0.8 else 'High' if normalized_score >= 0.6 else 'Moderate' if normalized_score >= 0.4 else 'Low' if normalized_score >= 0.2 else 'Poor',
            'details': details
        }

    def run_analysis(self, genome_file: str = "Genome.txt") -> Dict:
        """Run comprehensive longevity biomarkers analysis"""
        print("🧬 LONGEVITY BIOMARKERS ANALYSIS")
        print("=" * 60)
        print("Focusing on POSITIVE longevity-promoting markers")
        print("=" * 60)
        
        # Load genome data
        genome_data = self.load_genome_data(genome_file)
        
        # Run all analyses
        longevity_results = self.analyze_longevity_biomarkers(genome_data)
        centenarian_results = self.analyze_centenarian_traits(genome_data)
        anti_aging_results = self.analyze_anti_aging_markers(genome_data)
        probability_results = self.calculate_longevity_probability(genome_data)
        
        # Compile results
        results = {
            'longevity': longevity_results,
            'centenarian': centenarian_results,
            'anti_aging': anti_aging_results,
            'probability': probability_results
        }
        
        # Print detailed results
        print(f"\n📊 LONGEVITY BIOMARKERS ANALYSIS")
        print(f"Markers found: {longevity_results['markers_found']}/{longevity_results['total_markers']}")
        print(f"Protective alleles found: {longevity_results['protective_alleles_found']}")
        print(f"Longevity Score: {longevity_results['longevity_score']:.3f}")
        print(f"Longevity Level: {longevity_results['longevity_level']}")
        
        print(f"\n📊 CENTENARIAN TRAITS ANALYSIS")
        print(f"Markers found: {centenarian_results['markers_found']}/{centenarian_results['total_markers']}")
        print(f"Protective alleles found: {centenarian_results['protective_alleles_found']}")
        print(f"Centenarian Score: {centenarian_results['centenarian_score']:.3f}")
        print(f"Centenarian Potential: {centenarian_results['centenarian_potential']}")
        
        print(f"\n📊 ANTI-AGING MARKERS ANALYSIS")
        print(f"Markers found: {anti_aging_results['markers_found']}/{anti_aging_results['total_markers']}")
        print(f"Protective alleles found: {anti_aging_results['protective_alleles_found']}")
        print(f"Anti-Aging Score: {anti_aging_results['anti_aging_score']:.3f}")
        print(f"Anti-Aging Potential: {anti_aging_results['anti_aging_potential']}")
        
        print(f"\n🎯 POSITIVE LONGEVITY PROBABILITY SUMMARY")
        print(f"Total markers analyzed: {probability_results['total_markers']}")
        print(f"Markers used: {probability_results['markers_used']}")
        print(f"Protective alleles found: {probability_results['protective_alleles_found']}")
        print(f"Longevity Score: {probability_results['longevity_score']:.3f}")
        print(f"Overall Longevity Potential: {probability_results['overall_longevity_potential']}")
        
        print(f"\n📈 POSITIVE PROBABILITY REPORT")
        print(f"Longevity Probability: {probability_results['longevity_probability']:.1%}")
        print(f"Centenarian Probability: {probability_results['centenarian_probability']:.1%}")
        print(f"Anti-Aging Probability: {probability_results['anti_aging_probability']:.1%}")
        print(f"Healthy Aging Probability: {probability_results['healthy_aging_probability']:.1%}")
        
        return results

def main():
    """Main function to run the analysis"""
    analyzer = LongevityBiomarkersAnalyzer()
    results = analyzer.run_analysis()
    
    print(f"\n✅ Analysis complete!")
    print(f"📋 Summary: {results['probability']['overall_longevity_potential']} longevity potential")
    print(f"🎯 Focus on positive biomarkers for longevity")

if __name__ == "__main__":
    main() 