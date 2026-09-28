if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Carrier screening analysis for inherited conditions."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data
from config.carrier_screening_snps import CARRIER_SCREENING_SNPS, POPULATION_CARRIER_FREQUENCIES
from typing import Dict, Any, Optional


def analyze_carrier_status(genome_data) -> Dict[str, Any]:
    """
    Analyze carrier status for inherited conditions from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing carrier screening results
    """
    results = {
        'conditions': {},
        'carriers_detected': [],
        'summary': {
            'total_markers': len(CARRIER_SCREENING_SNPS),
            'markers_found': 0,
            'carriers_found': 0
        }
    }
    
    # Group by condition
    conditions = {}
    for rsid, data in CARRIER_SCREENING_SNPS.items():
        condition = data['condition']
        if condition not in conditions:
            conditions[condition] = []
        conditions[condition].append((rsid, data))
    
    # Analyze each condition
    for condition, variants in conditions.items():
        condition_results = {
            'condition': condition,
            'variants': [],
            'carrier_status': 'Not a carrier',
            'risk_level': 'Low'
        }
        
        for rsid, variant_data in variants:
            genotype = genome_data.get_genotype(rsid)
            
            variant_result = {
                'rsid': rsid,
                'gene': variant_data['gene'],
                'variant': variant_data['variant'],
                'found': genotype is not None,
            }
            
            if genotype:
                results['summary']['markers_found'] += 1
                risk_allele = variant_data['risk_allele']
                
                # Check for variant
                allele_count = genotype.count(risk_allele)
                is_homozygous = allele_count == 2
                is_heterozygous = allele_count == 1
                
                variant_result['genotype'] = genotype
                variant_result['risk_allele'] = risk_allele
                variant_result['allele_count'] = allele_count
                variant_result['is_carrier'] = is_heterozygous
                variant_result['is_affected'] = is_homozygous
                
                if is_homozygous:
                    condition_results['carrier_status'] = 'Affected (homozygous)'
                    condition_results['risk_level'] = 'High'
                    variant_result['status'] = 'Affected'
                elif is_heterozygous:
                    condition_results['carrier_status'] = 'Carrier (heterozygous)'
                    condition_results['risk_level'] = 'Carrier'
                    variant_result['status'] = 'Carrier'
                    results['summary']['carriers_found'] += 1
                    if condition not in results['carriers_detected']:
                        results['carriers_detected'].append(condition)
                else:
                    variant_result['status'] = 'Not a carrier'
                
                variant_result['description'] = variant_data['description']
                variant_result['prevalence'] = variant_data.get('prevalence', 'Unknown')
            else:
                variant_result['genotype'] = None
                variant_result['status'] = 'Not tested'
            
            condition_results['variants'].append(variant_result)
        
        results['conditions'][condition] = condition_results
    
    return results


def print_carrier_screening_report(results: Dict[str, Any]):
    """Print formatted carrier screening report to console."""
    print("\n" + "="*60)
    print("CARRIER SCREENING ANALYSIS")
    print("="*60 + "\n")
    print("Inherited Condition Carrier Status")
    print("-"*60 + "\n")
    
    # Print each condition
    for condition, condition_data in results['conditions'].items():
        print("="*60)
        print(f"CONDITION: {condition}")
        print("="*60)
        print(f"Carrier Status: {condition_data['carrier_status']}")
        print(f"Risk Level: {condition_data['risk_level']}")
        print()
        
        for variant in condition_data['variants']:
            print(f"Gene: {variant['gene']} ({variant['rsid']})")
            print(f"  Variant: {variant['variant']}")
            if variant['found']:
                print(f"  Genotype: {variant['genotype']}")
                print(f"  Status: {variant['status']}")
                if variant.get('allele_count') is not None:
                    print(f"  Allele Count: {variant['allele_count']}")
                print(f"  Description: {variant['description']}")
                if variant.get('prevalence') != 'Unknown':
                    print(f"  Population Prevalence: {variant['prevalence']}")
            else:
                print(f"  Status: NOT FOUND IN DATA")
            print()
    
    # Print carriers detected
    if results['carriers_detected']:
        print("="*60)
        print("CARRIERS DETECTED")
        print("="*60)
        for condition in results['carriers_detected']:
            print(f"- {condition}")
        print()
    
    # Summary
    summary = results['summary']
    print("="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total markers analyzed: {summary['total_markers']}")
    print(f"Markers found: {summary['markers_found']}")
    print(f"Carriers detected: {summary['carriers_found']}")
    
    # Important information
    print("\n" + "="*60)
    print("IMPORTANT INFORMATION")
    print("="*60)
    print("- Carrier status means you have one copy of a disease-causing variant")
    print("- Carriers are typically healthy but can pass the variant to children")
    print("- If both parents are carriers, child has 25% chance of being affected")
    print("- Consult genetic counselors for family planning decisions")
    print("- This screening may not detect all variants")
    print("- This is for educational purposes only, not medical advice")


def main(genome_file: str = "Genome.txt", output_json: Optional[str] = None):
    """
    Main function to run carrier screening analysis.
    
    Args:
        genome_file: Path to genome file
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze carrier status
        results = analyze_carrier_status(genome_data)
        
        # Print report
        print_carrier_screening_report(results)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'carrier_screening': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Carrier Screening Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, output_json=args.json)

