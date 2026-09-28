if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Sleep and chronotype analysis."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data
from config.sleep_snps import SLEEP_SNPS, CHRONOTYPE_INTERPRETATION
from typing import Dict, Any, Optional


def analyze_sleep(genome_data) -> Dict[str, Any]:
    """
    Analyze sleep and chronotype markers from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing sleep analysis results
    """
    results = {
        'sleep_traits': {},
        'chronotype_score': {'morning': 0, 'evening': 0},
        'summary': {
            'total_markers': len(SLEEP_SNPS),
            'markers_found': 0,
            'variants_detected': 0
        }
    }
    
    # Analyze each sleep trait
    for rsid, snp_data in SLEEP_SNPS.items():
        trait = snp_data['trait']
        
        if trait not in results['sleep_traits']:
            results['sleep_traits'][trait] = {
                'trait': trait,
                'variants': []
            }
        
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
            variant_result['interpretation'] = snp_data.get('interpretation', '')
            
            if has_variant:
                results['summary']['variants_detected'] += 1
                
                # Score chronotype
                if 'morning' in snp_data['description'].lower() or 'morning' in snp_data.get('interpretation', '').lower():
                    results['chronotype_score']['morning'] += 1
                elif 'evening' in snp_data['description'].lower() or 'evening' in snp_data.get('interpretation', '').lower() or 'night owl' in snp_data.get('interpretation', '').lower():
                    results['chronotype_score']['evening'] += 1
        else:
            variant_result['genotype'] = None
            variant_result['variant_present'] = False
        
        results['sleep_traits'][trait]['variants'].append(variant_result)
    
    # Determine chronotype
    morning_score = results['chronotype_score']['morning']
    evening_score = results['chronotype_score']['evening']
    
    if morning_score > evening_score:
        results['predicted_chronotype'] = 'Morning'
    elif evening_score > morning_score:
        results['predicted_chronotype'] = 'Evening'
    else:
        results['predicted_chronotype'] = 'Intermediate'
    
    return results


def print_sleep_report(results: Dict[str, Any]):
    """Print formatted sleep analysis report to console."""
    print("\n" + "="*60)
    print("SLEEP & CHRONOTYPE ANALYSIS")
    print("="*60 + "\n")
    
    # Print chronotype prediction
    print("="*60)
    print("PREDICTED CHRONOTYPE")
    print("="*60)
    chronotype = results['predicted_chronotype']
    print(f"Chronotype: {chronotype}")
    print(f"Interpretation: {CHRONOTYPE_INTERPRETATION.get(chronotype, 'Flexible sleep schedule')}")
    print(f"Morning Score: {results['chronotype_score']['morning']}")
    print(f"Evening Score: {results['chronotype_score']['evening']}")
    print()
    
    # Print trait analysis
    for trait, trait_data in results['sleep_traits'].items():
        print("="*60)
        print(f"TRAIT: {trait}")
        print("="*60)
        
        for variant in trait_data['variants']:
            print(f"Gene: {variant['gene']} ({variant['rsid']})")
            if variant['found']:
                print(f"  Genotype: {variant['genotype']}")
                print(f"  Variant Present: {variant['variant_present']}")
                print(f"  Description: {variant['description']}")
                if variant.get('interpretation'):
                    print(f"  Interpretation: {variant['interpretation']}")
            else:
                print(f"  Status: NOT FOUND IN DATA")
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
    print("SLEEP RECOMMENDATIONS")
    print("="*60)
    print("- Maintain consistent sleep schedule aligned with your chronotype")
    print("- Morning types: Schedule important tasks in morning")
    print("- Evening types: Allow flexibility for later sleep times when possible")
    print("- Practice good sleep hygiene regardless of chronotype")
    print("- Consult sleep specialists if experiencing sleep disorders")
    print("\n- This is for educational purposes only, not medical advice")


def main(genome_file: str = "Genome.txt", output_json: Optional[str] = None):
    """
    Main function to run sleep analysis.
    
    Args:
        genome_file: Path to genome file
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze sleep
        results = analyze_sleep(genome_data)
        
        # Print report
        print_sleep_report(results)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'sleep': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Sleep & Chronotype Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, output_json=args.json)

