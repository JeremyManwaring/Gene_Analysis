if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Comprehensive disease risk analysis using shared utilities and centralized SNP database."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data, get_allele_count, calculate_prs
from config import DISEASE_SNPS, DISEASE_STATS, DISEASE_PRS_MODELS
from typing import Dict, Any, Optional


def analyze_disease_risk(genome_data, detailed: bool = True) -> Dict[str, Any]:
    """
    Analyze disease risk from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
        detailed: If True, include detailed statistics and PRS analysis
    
    Returns:
        Dictionary containing disease risk analysis results
    """
    results = {
        'disease_risks': [],
        'summary': {
            'total_snps': len(DISEASE_SNPS),
            'snps_found': 0,
            'risk_alleles_found': 0,
        }
    }
    
    # Group SNPs by population
    populations = {}
    for rsid, data in DISEASE_SNPS.items():
        pop = data['population']
        if pop not in populations:
            populations[pop] = []
        populations[pop].append((rsid, data))
    
    # Analyze each population group
    for population in ['Global', 'East Asian']:
        if population not in populations:
            continue
        
        population_results = {
            'population': population,
            'diseases': []
        }
        
        for rsid, data in populations[population]:
            genotype = genome_data.get_genotype(rsid)
            
            disease_result = {
                'rsid': rsid,
                'disease': data['disease'],
                'gene': data.get('gene', 'Unknown'),
                'description': data['description'],
                'found': genotype is not None,
            }
            
            if genotype:
                results['summary']['snps_found'] += 1
                has_risk = any(allele in genotype for allele in data['risk'])
                
                disease_result['genotype'] = genotype
                disease_result['risk_alleles'] = data['risk']
                disease_result['risk_detected'] = has_risk
                
                if has_risk:
                    results['summary']['risk_alleles_found'] += 1
                
                # Add statistics if available
                if detailed and rsid in DISEASE_STATS:
                    stats = DISEASE_STATS[rsid]
                    disease_result['statistics'] = {
                        'population_frequency': stats['population_frequency'],
                        'odds_ratio': stats['odds_ratio'],
                        'confidence_interval': stats['confidence_interval'],
                    }
            else:
                disease_result['genotype'] = None
                disease_result['risk_detected'] = None
            
            population_results['diseases'].append(disease_result)
        
        results['disease_risks'].append(population_results)
    
    # Calculate coverage
    if results['summary']['total_snps'] > 0:
        results['summary']['coverage'] = (
            results['summary']['snps_found'] / results['summary']['total_snps']
        ) * 100
    else:
        results['summary']['coverage'] = 0.0
    
    # Add PRS analysis if detailed
    if detailed:
        results['prs_scores'] = analyze_prs(genome_data)
    
    return results


def analyze_prs(genome_data) -> Dict[str, Any]:
    """
    Calculate polygenic risk scores for various diseases.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing PRS scores for each disease model
    """
    prs_results = {}
    
    for trait, snps in DISEASE_PRS_MODELS.items():
        # Prepare genotype data and weights
        genotype_dict = {}
        snp_weights = {}
        
        for rsid, (risk_allele, weight, source) in snps.items():
            genotype = genome_data.get_genotype(rsid)
            if genotype:
                genotype_dict[rsid] = genotype
                snp_weights[rsid] = (risk_allele, weight)
        
        # Calculate PRS
        prs, snps_found = calculate_prs(genotype_dict, snp_weights, use_effect=False)
        
        prs_results[trait] = {
            'prs_score': round(prs, 2),
            'snps_found': snps_found,
            'total_snps': len(snps),
            'coverage': (snps_found / len(snps) * 100) if snps else 0.0,
        }
    
    return prs_results


def print_disease_report(results: Dict[str, Any], detailed: bool = True):
    """Print formatted disease risk report to console."""
    print("\n" + "="*60)
    print("COMPREHENSIVE DISEASE RISK ANALYSIS")
    print("="*60 + "\n")
    
    # Print disease risks by population
    for pop_result in results['disease_risks']:
        print("="*60)
        print(f"POPULATION: {pop_result['population']}")
        print("="*60 + "\n")
        
        for disease in pop_result['diseases']:
            print("-"*50)
            print(f"Disease: {disease['disease']}")
            print("-"*50)
            print(f"SNP: {disease['rsid']}")
            print(f"Gene: {disease['gene']}")
            
            if disease['found']:
                print(f"Your Genotype: {disease['genotype']}")
                print(f"Risk Allele(s): {', '.join(disease['risk_alleles'])}")
                status = "RISK DETECTED" if disease['risk_detected'] else "NO RISK"
                print(f"Status: {status}")
            else:
                print("SNP: NOT FOUND IN DATA")
            
            print(f"Description: {disease['description']}")
            
            if detailed and 'statistics' in disease:
                stats = disease['statistics']
                print(f"Population Frequency: {stats['population_frequency']}")
                print(f"Odds Ratio: {stats['odds_ratio']}")
                print(f"Confidence Interval: {stats['confidence_interval']}")
            
            print()
    
    # Print summary
    summary = results['summary']
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Total SNPs analyzed: {summary['total_snps']}")
    print(f"SNPs found in your data: {summary['snps_found']}")
    print(f"Risk alleles found: {summary['risk_alleles_found']}")
    print(f"Data coverage: {summary['coverage']:.1f}%")
    
    # Print PRS scores if available
    if 'prs_scores' in results:
        print("\n" + "="*60)
        print("POLYGENIC RISK SCORES (PRS)")
        print("="*60 + "\n")
        
        for trait, prs_data in results['prs_scores'].items():
            print("="*50)
            print(f"TRAIT: {trait}")
            print("="*50)
            print(f"PRS Score: {prs_data['prs_score']}")
            print(f"Based on {prs_data['snps_found']}/{prs_data['total_snps']} SNPs found")
            print(f"Coverage: {prs_data['coverage']:.1f}%")
            print()
    
    # Print interpretation
    print("="*60)
    print("INTERPRETATION")
    print("="*60)
    print("- RISK DETECTED: You carry risk alleles for this condition")
    print("- NO RISK: You do not carry the known risk alleles")
    print("- NOT FOUND: This SNP was not present in your genetic data")
    if 'prs_scores' in results:
        print("- PRS: Higher scores = higher genetic risk compared to population average")
    print("- These results are for educational purposes only")
    print("- Consult healthcare professionals for medical advice")


def main(genome_file: str = "Genome.txt", detailed: bool = True, output_json: Optional[str] = None):
    """
    Main function to run disease analysis.
    
    Args:
        genome_file: Path to genome file
        detailed: Include detailed statistics and PRS
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze disease risk
        results = analyze_disease_risk(genome_data, detailed=detailed)
        
        # Print report
        print_disease_report(results, detailed=detailed)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'disease': results}, filename=output_json)
            print(f"\nResults exported to: {output_json}")
        
        return results
        
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error during analysis: {e}")
        sys.exit(1)


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Disease Risk Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--simple", action="store_true", help="Simple output (no detailed stats)")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, detailed=not args.simple, output_json=args.json)

