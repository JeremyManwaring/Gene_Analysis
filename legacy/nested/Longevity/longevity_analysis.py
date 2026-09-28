if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Comprehensive longevity analysis with probability estimates using shared utilities."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data, genotype_effect, normalize_score
from config import LONGEVITY_MARKERS, LONGEVITY_TOTAL_MARKERS
from typing import Dict, Any, Optional


def analyze_longevity(genome_data, comprehensive: bool = True) -> Dict[str, Any]:
    """
    Analyze longevity markers from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
        comprehensive: If True, use comprehensive marker set
    
    Returns:
        Dictionary containing longevity analysis results
    """
    results = {
        'markers': [],
        'summary': {
            'total_markers': len(LONGEVITY_MARKERS),
            'markers_found': 0,
        }
    }
    
    raw_score = 0.0
    
    # Analyze each marker
    for rsid, info in LONGEVITY_MARKERS.items():
        genotype = genome_data.get_genotype(rsid)
        
        marker_result = {
            'rsid': rsid,
            'description': info['description'],
            'found': genotype is not None,
        }
        
        if genotype:
            results['summary']['markers_found'] += 1
            effect = genotype_effect(genotype, info['risk'][0])
            contribution = effect * info['weight']
            raw_score += contribution
            
            marker_result['genotype'] = genotype
            marker_result['risk_allele'] = info['risk'][0]
            marker_result['weight'] = info['weight']
            marker_result['effect'] = effect
            marker_result['contribution'] = round(contribution, 3)
        else:
            marker_result['genotype'] = None
            marker_result['contribution'] = 0.0
        
        results['markers'].append(marker_result)
    
    # Calculate normalized score
    normalized_score = normalize_score(raw_score, LONGEVITY_TOTAL_MARKERS)
    aging_rate_probability = 1 - normalized_score
    longevity_probability = normalized_score
    
    # Determine risk level
    if normalized_score < 0.4:
        risk_level = 'High'
    elif normalized_score < 0.7:
        risk_level = 'Moderate'
    else:
        risk_level = 'Low'
    
    results['summary'].update({
        'raw_score': round(raw_score, 3),
        'normalized_score': round(normalized_score, 3),
        'aging_rate_probability': round(aging_rate_probability, 3),
        'longevity_probability': round(longevity_probability, 3),
        'cellular_health_probability': round(normalized_score, 3),
        'risk_level': risk_level,
        'coverage': (results['summary']['markers_found'] / len(LONGEVITY_MARKERS) * 100) if LONGEVITY_MARKERS else 0.0,
    })
    
    return results


def print_longevity_report(results: Dict[str, Any], comprehensive: bool = True):
    """Print formatted longevity report to console."""
    print("\n" + "="*60)
    if comprehensive:
        print("COMPREHENSIVE LONGEVITY REPORT")
    else:
        print("LONGEVITY MARKER SUMMARY")
    print("="*60 + "\n")
    
    # Print marker details
    for marker in results['markers']:
        if marker['found']:
            print(f"SNP: {marker['rsid']:<12} Genotype: {marker['genotype']:<6} "
                  f"Effect: {marker['contribution']:.3f}  {marker['description']}")
        else:
            print(f"SNP: {marker['rsid']:<12} NOT FOUND  {marker['description']}")
    
    # Print summary
    summary = results['summary']
    print(f"\nMarkers found: {summary['markers_found']}/{summary['total_markers']}")
    print(f"Polygenic score: {summary['normalized_score']:.3f}")
    print(f"Aging rate probability: {summary['aging_rate_probability']:.3f} ({summary['aging_rate_probability']*100:.1f}%)")
    print(f"Longevity probability: {summary['longevity_probability']:.3f} ({summary['longevity_probability']*100:.1f}%)")
    print(f"Cellular health probability: {summary['cellular_health_probability']:.3f} ({summary['cellular_health_probability']*100:.1f}%)")
    print(f"Risk level: {summary['risk_level']}")
    print(f"Data coverage: {summary['coverage']:.1f}%")
    
    # Print interpretation
    print("\n" + "="*60)
    print("INTERPRETATION")
    print("="*60)
    print(f"• Higher polygenic score = better aging genetics")
    print(f"• Your score: {summary['normalized_score']:.3f}")
    print(f"• Risk level: {summary['risk_level']}")
    if summary['risk_level'] == 'High':
        print(f"  - High risk of accelerated aging ({summary['aging_rate_probability']*100:.1f}%)")
        print(f"  - Low probability of extended lifespan ({summary['longevity_probability']*100:.1f}%)")
        print(f"  - Low cellular health maintenance ({summary['cellular_health_probability']*100:.1f}%)")
    elif summary['risk_level'] == 'Moderate':
        print(f"  - Moderate aging risk")
        print(f"  - Moderate longevity potential")
    else:
        print(f"  - Lower aging risk")
        print(f"  - Higher longevity potential")
    print("\n• These results are for educational purposes only")
    print("• Consult healthcare professionals for medical advice")


def main(genome_file: str = "Genome.txt", comprehensive: bool = True, output_json: Optional[str] = None):
    """
    Main function to run longevity analysis.
    
    Args:
        genome_file: Path to genome file
        comprehensive: Use comprehensive marker set
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze longevity
        results = analyze_longevity(genome_data, comprehensive=comprehensive)
        
        # Print report
        print_longevity_report(results, comprehensive=comprehensive)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'longevity': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Longevity Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--simple", action="store_true", help="Simple analysis (fewer markers)")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, comprehensive=not args.simple, output_json=args.json)

