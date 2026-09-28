if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Comprehensive ancestry analysis using shared utilities."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data
from config import ANCESTRY_CATEGORIES
from typing import Dict, Any, Optional


def analyze_ancestry(genome_data) -> Dict[str, Any]:
    """
    Analyze ancestry-related markers from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing ancestry analysis results
    """
    results = {
        'categories': {},
        'summary': {},
        'population_markers': {}
    }
    
    all_results = {}
    category_counts = {}
    
    # Analyze each category
    for category, snps in ANCESTRY_CATEGORIES.items():
        category_results = {
            'category': category,
            'snps': []
        }
        
        found_count = 0
        
        for rsid, description in snps.items():
            genotype = genome_data.get_genotype(rsid)
            
            snp_result = {
                'rsid': rsid,
                'description': description,
                'found': genotype is not None,
            }
            
            if genotype:
                found_count += 1
                snp_result['genotype'] = genotype
                all_results[rsid] = (genotype, description, category)
            else:
                snp_result['genotype'] = None
            
            category_results['snps'].append(snp_result)
        
        category_counts[category] = (found_count, len(snps))
        results['categories'][category] = category_results
    
    # Calculate summary statistics
    total_found = sum(count for count, _ in category_counts.values())
    total_snps = sum(total for _, total in category_counts.values())
    
    results['summary'] = {
        'total_snps': total_snps,
        'snps_found': total_found,
        'snps_missing': total_snps - total_found,
        'coverage': (total_found / total_snps * 100) if total_snps > 0 else 0.0,
        'category_breakdown': {
            category: {
                'found': found,
                'total': total,
                'coverage': (found / total * 100) if total > 0 else 0.0
            }
            for category, (found, total) in category_counts.items()
        }
    }
    
    # Population-specific analysis
    east_asian_markers = ['rs3827760']
    european_markers = [
        'rs1426654', 'rs16891982', 'rs12913832', 'rs1805007',
        'rs1805008', 'rs885479', 'rs12203592', 'rs6058017'
    ]
    african_markers = ['rs1426654', 'rs16891982', 'rs12913832']
    
    east_asian_count = sum(1 for rsid in all_results if rsid in east_asian_markers)
    european_count = sum(1 for rsid in all_results if rsid in european_markers)
    african_count = sum(1 for rsid in all_results if rsid in african_markers)
    
    results['population_markers'] = {
        'east_asian': {
            'markers_found': east_asian_count,
            'total_markers': len(east_asian_markers)
        },
        'european': {
            'markers_found': european_count,
            'total_markers': len(european_markers)
        },
        'african': {
            'markers_found': african_count,
            'total_markers': len(african_markers)
        }
    }
    
    # Key findings
    key_findings = []
    
    if 'rs3827760' in all_results:
        genotype = all_results['rs3827760'][0]
        if genotype in ('AA', 'AG'):
            key_findings.append("✓ East Asian hair thickness variant detected")
        else:
            key_findings.append("✗ East Asian hair thickness variant not detected")
    
    if 'rs671' in all_results:
        genotype = all_results['rs671'][0]
        if genotype in ('AA', 'AG'):
            key_findings.append("✓ Asian flush (ALDH2 deficiency) variant detected")
        else:
            key_findings.append("✗ Asian flush variant not detected")
    
    if 'rs4988235' in all_results:
        genotype = all_results['rs4988235'][0]
        if 'C' in genotype:
            key_findings.append("✓ Lactose intolerance variant detected")
        else:
            key_findings.append("✗ Lactose intolerance variant not detected")
    
    results['key_findings'] = key_findings
    
    return results


def print_ancestry_report(results: Dict[str, Any]):
    """Print formatted ancestry analysis report to console."""
    print("\n" + "="*60)
    print("COMPREHENSIVE ANCESTRY SUMMARY")
    print("="*60 + "\n")
    
    # Print each category
    for category, category_data in results['categories'].items():
        print("="*60)
        print(f"CATEGORY: {category}")
        print("="*60)
        
        for snp in category_data['snps']:
            print("-"*50)
            print(f"SNP: {snp['rsid']}")
            if snp['found']:
                print(f"Genotype: {snp['genotype']}")
            else:
                print("SNP: NOT FOUND IN DATA")
            print(f"Function: {snp['description']}")
            print()
    
    # Summary statistics
    summary = results['summary']
    print("="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Total SNPs analyzed: {summary['total_snps']}")
    print(f"SNPs found in your data: {summary['snps_found']}")
    print(f"SNPs missing: {summary['snps_missing']}")
    print(f"Data coverage: {summary['coverage']:.1f}%")
    
    # Category breakdown
    print("\n" + "="*60)
    print("CATEGORY BREAKDOWN")
    print("="*60)
    
    for category, breakdown in summary['category_breakdown'].items():
        print(f"{category:<25} {breakdown['found']}/{breakdown['total']} SNPs "
              f"({breakdown['coverage']:.1f}% coverage)")
    
    # Population-specific analysis
    print("\n" + "="*60)
    print("POPULATION-SPECIFIC ANALYSIS")
    print("="*60)
    
    pop_markers = results['population_markers']
    print(f"East Asian markers found: {pop_markers['east_asian']['markers_found']}")
    print(f"European markers found: {pop_markers['european']['markers_found']}")
    print(f"African markers found: {pop_markers['african']['markers_found']}")
    
    # Key findings
    if results['key_findings']:
        print("\n" + "="*60)
        print("KEY FINDINGS")
        print("="*60)
        for finding in results['key_findings']:
            print(finding)
    
    # Interpretation
    print("\n" + "="*60)
    print("INTERPRETATION")
    print("="*60)
    print("- This analysis examines genetic markers that differ between populations")
    print("- Higher counts in a category suggest genetic similarity to that population")
    print("- Remember: Genetic ancestry ≠ Cultural identity")
    print("- These results are statistical associations, not definitive ancestry")
    print("- Individual genetic variation is complex and doesn't fit simple categories")
    print("- Your cultural identity is valid regardless of genetic markers")


def main(genome_file: str = "Genome.txt", output_json: Optional[str] = None):
    """
    Main function to run ancestry analysis.
    
    Args:
        genome_file: Path to genome file
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze ancestry
        results = analyze_ancestry(genome_data)
        
        # Print report
        print_ancestry_report(results)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'ancestry': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Ancestry Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, output_json=args.json)

