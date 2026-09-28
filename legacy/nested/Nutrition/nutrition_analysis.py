if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Nutrition and vitamin metabolism analysis."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data
from config.nutrition_snps import NUTRITION_SNPS, VITAMIN_RECOMMENDATIONS
from typing import Dict, Any, Optional


def analyze_nutrition(genome_data) -> Dict[str, Any]:
    """
    Analyze nutrition and vitamin metabolism markers from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing nutrition analysis results
    """
    results = {
        'vitamin_traits': {},
        'recommendations': {},
        'summary': {
            'total_markers': len(NUTRITION_SNPS),
            'markers_found': 0,
            'variants_detected': 0
        }
    }
    
    # Group by trait/vitamin
    traits = {}
    for rsid, data in NUTRITION_SNPS.items():
        trait = data['trait']
        if trait not in traits:
            traits[trait] = []
        traits[trait].append((rsid, data))
    
    # Analyze each trait
    for trait, snps in traits.items():
        trait_results = {
            'trait': trait,
            'variants': [],
            'recommendations': []
        }
        
        for rsid, snp_data in snps:
            genotype = genome_data.get_genotype(rsid)
            
            variant_result = {
                'rsid': rsid,
                'gene': snp_data['gene'],
                'found': genotype is not None,
            }
            
            if genotype:
                results['summary']['markers_found'] += 1
                risk_allele = snp_data['risk_allele']
                has_variant = risk_allele in genotype
                
                variant_result['genotype'] = genotype
                variant_result['risk_allele'] = risk_allele
                variant_result['variant_present'] = has_variant
                variant_result['description'] = snp_data['description']
                
                if has_variant:
                    results['summary']['variants_detected'] += 1
                    variant_result['recommendation'] = snp_data.get('recommendation', '')
                    if variant_result['recommendation']:
                        trait_results['recommendations'].append(variant_result['recommendation'])
            else:
                variant_result['genotype'] = None
                variant_result['variant_present'] = False
            
            trait_results['variants'].append(variant_result)
        
        results['vitamin_traits'][trait] = trait_results
    
    # Add general vitamin recommendations
    for vitamin, rec_data in VITAMIN_RECOMMENDATIONS.items():
        relevant_snps = rec_data['snps']
        found_variants = []
        
        for rsid in relevant_snps:
            if rsid in NUTRITION_SNPS:
                genotype = genome_data.get_genotype(rsid)
                if genotype:
                    snp_data = NUTRITION_SNPS[rsid]
                    if snp_data['risk_allele'] in genotype:
                        found_variants.append(rsid)
        
        if found_variants:
            results['recommendations'][vitamin] = {
                'variants_found': found_variants,
                'recommendation': rec_data['general']
            }
    
    return results


def print_nutrition_report(results: Dict[str, Any]):
    """Print formatted nutrition analysis report to console."""
    print("\n" + "="*60)
    print("NUTRITION & VITAMIN METABOLISM ANALYSIS")
    print("="*60 + "\n")
    
    # Print trait analysis
    for trait, trait_data in results['vitamin_traits'].items():
        print("="*60)
        print(f"TRAIT: {trait}")
        print("="*60)
        
        for variant in trait_data['variants']:
            print(f"Gene: {variant['gene']} ({variant['rsid']})")
            if variant['found']:
                print(f"  Genotype: {variant['genotype']}")
                print(f"  Variant Present: {variant['variant_present']}")
                print(f"  Description: {variant['description']}")
                if variant.get('recommendation'):
                    print(f"  Recommendation: {variant['recommendation']}")
            else:
                print(f"  Status: NOT FOUND IN DATA")
            print()
        
        if trait_data['recommendations']:
            print("Recommendations:")
            for rec in trait_data['recommendations']:
                print(f"  - {rec}")
            print()
    
    # Print vitamin recommendations
    if results['recommendations']:
        print("="*60)
        print("VITAMIN RECOMMENDATIONS")
        print("="*60 + "\n")
        
        for vitamin, rec_data in results['recommendations'].items():
            print(f"{vitamin}:")
            print(f"  Variants Found: {', '.join(rec_data['variants_found'])}")
            print(f"  Recommendation: {rec_data['recommendation']}")
            print()
    
    # Summary
    summary = results['summary']
    print("="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total markers analyzed: {summary['total_markers']}")
    print(f"Markers found: {summary['markers_found']}")
    print(f"Variants detected: {summary['variants_detected']}")
    
    print("\n" + "="*60)
    print("IMPORTANT NOTES")
    print("="*60)
    print("- These are genetic predispositions, not diagnoses")
    print("- Consult healthcare providers before taking supplements")
    print("- Blood tests may be needed to confirm deficiencies")
    print("- Diet and lifestyle also affect nutrient levels")
    print("- This is for educational purposes only, not medical advice")


def main(genome_file: str = "Genome.txt", output_json: Optional[str] = None):
    """
    Main function to run nutrition analysis.
    
    Args:
        genome_file: Path to genome file
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze nutrition
        results = analyze_nutrition(genome_data)
        
        # Print report
        print_nutrition_report(results)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'nutrition': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Nutrition & Vitamin Metabolism Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, output_json=args.json)

