if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Comprehensive fitness and athletic performance analysis using shared utilities."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data, get_allele_count, calculate_prs, calculate_probability
from utils.output_formatter import format_snp_table
from config import FITNESS_SNPS
from typing import Dict, Any, Optional, List, Tuple


def analyze_fitness(genome_data) -> Dict[str, Any]:
    """
    Analyze fitness and athletic traits from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing fitness analysis results
    """
    results = {
        'traits': {},
        'summary': {}
    }
    
    trait_probabilities = {}
    
    for trait, snps in FITNESS_SNPS.items():
        # Prepare genotype data and weights
        genotype_dict = {}
        snp_weights = {}
        snp_details = []
        
        for rsid, data in snps.items():
            genotype = genome_data.get_genotype(rsid)
            if genotype:
                genotype_dict[rsid] = genotype
                snp_weights[rsid] = (data['risk'][0], data['weight'])
                count = get_allele_count(genotype, data['risk'][0])
                snp_details.append((rsid, genotype, data['risk'][0], data['weight'], count, data['description']))
            else:
                snp_details.append((rsid, 'Not found', data['risk'][0], data['weight'], 0, data['description']))
        
        # Calculate PRS
        prs, snps_found = calculate_prs(genotype_dict, snp_weights, use_effect=False)
        
        # Calculate probability
        probability = calculate_probability(prs)
        trait_probabilities[trait] = probability
        
        # Determine interpretation
        if probability >= 70:
            interpretation = "HIGH"
        elif probability >= 50:
            interpretation = "MODERATE"
        else:
            interpretation = "LOW"
        
        results['traits'][trait] = {
            'prs_score': round(prs, 2),
            'probability': probability,
            'interpretation': interpretation,
            'snps_found': snps_found,
            'total_snps': len(snps),
            'coverage': (snps_found / len(snps) * 100) if snps else 0.0,
            'snp_details': snp_details,
        }
    
    results['summary'] = {
        'total_traits': len(FITNESS_SNPS),
        'trait_probabilities': trait_probabilities,
    }
    
    return results


def get_trait_interpretation(trait: str, probability: float) -> str:
    """Get detailed interpretation for a fitness trait."""
    if trait == 'Muscle Fiber Type':
        if probability >= 70:
            return (
                f"{probability}% probability of FAST-TWITCH muscle dominance\n"
                "  → Better suited for power sports (sprinting, weightlifting)"
            )
        elif probability >= 50:
            return (
                f"{probability}% probability of MIXED muscle fiber type\n"
                "  → Balanced for both power and endurance"
            )
        else:
            return (
                f"{probability}% probability of SLOW-TWITCH muscle dominance\n"
                "  → Better suited for endurance sports (distance running, cycling)"
            )
    
    elif trait == 'Endurance Capacity':
        if probability >= 70:
            return (
                f"{probability}% probability of HIGH endurance capacity\n"
                "  → Excellent for long-distance events"
            )
        elif probability >= 50:
            return (
                f"{probability}% probability of MODERATE endurance capacity\n"
                "  → Good baseline for endurance training"
            )
        else:
            return (
                f"{probability}% probability of LOWER endurance capacity\n"
                "  → May need more training for endurance events"
            )
    
    elif trait == 'Power & Strength':
        if probability >= 70:
            return (
                f"{probability}% probability of HIGH power potential\n"
                "  → Excellent for explosive movements"
            )
        elif probability >= 50:
            return (
                f"{probability}% probability of MODERATE power potential\n"
                "  → Good baseline for strength training"
            )
        else:
            return (
                f"{probability}% probability of LOWER power potential\n"
                "  → May need more focus on strength training"
            )
    
    elif trait == 'Injury Risk':
        if probability >= 70:
            return (
                f"{probability}% probability of HIGHER injury risk\n"
                "  → Focus on injury prevention and proper form"
            )
        elif probability >= 50:
            return (
                f"{probability}% probability of MODERATE injury risk\n"
                "  → Standard injury prevention recommended"
            )
        else:
            return (
                f"{probability}% probability of LOWER injury risk\n"
                "  → Good genetic resilience, but still practice safety"
            )
    
    elif trait == 'Recovery Rate':
        if probability >= 70:
            return (
                f"{probability}% probability of FAST recovery\n"
                "  → Can handle higher training volumes"
            )
        elif probability >= 50:
            return (
                f"{probability}% probability of MODERATE recovery\n"
                "  → Standard recovery protocols recommended"
            )
        else:
            return (
                f"{probability}% probability of SLOWER recovery\n"
                "  → May need more recovery time between sessions"
            )
    
    else:
        if probability >= 70:
            return f"{probability}% probability of HIGH {trait.lower()}"
        elif probability >= 50:
            return f"{probability}% probability of MODERATE {trait.lower()}"
        else:
            return f"{probability}% probability of LOWER {trait.lower()}"


def print_fitness_report(results: Dict[str, Any]):
    """Print formatted fitness analysis report to console."""
    print("\n" + "="*60)
    print("FITNESS & ATHLETIC PERFORMANCE ANALYSIS")
    print("="*60 + "\n")
    
    # Print each trait
    for trait, trait_data in results['traits'].items():
        print("="*60)
        print(f"TRAIT: {trait}")
        print("="*60)
        
        # Print SNP table
        headers = ["SNP", "Genotype", "Risk", "Weight", "Count", "Description"]
        print(format_snp_table(trait_data['snp_details'], headers=headers))
        print()
        
        print(f"PRS for {trait}: {trait_data['prs_score']:.2f} "
              f"(based on {trait_data['snps_found']}/{trait_data['total_snps']} SNPs found)")
        print(f"Probability Score: {trait_data['probability']}%")
        print()
    
    # Summary
    print("="*60)
    print("SUMMARY & PROBABILITY SCORES")
    print("="*60)
    
    print(f"\n{'Trait':<25}{'Probability':<15}{'Interpretation'}")
    print("-"*60)
    
    for trait, trait_data in results['traits'].items():
        print(f"{trait:<25}{trait_data['probability']:<15}%{trait_data['interpretation']}")
    
    # Detailed interpretation
    print("\n" + "="*60)
    print("DETAILED INTERPRETATION")
    print("="*60)
    
    for trait, trait_data in results['traits'].items():
        print(f"\n{trait}:")
        print(get_trait_interpretation(trait, trait_data['probability']))
    
    # Recommendations
    print("\n" + "="*60)
    print("TRAINING RECOMMENDATIONS")
    print("="*60)
    print("- These probabilities indicate genetic predispositions, not limitations")
    print("- Training can significantly improve performance regardless of genetics")
    print("- Focus on your strengths while developing areas of opportunity")
    print("- Consult with fitness professionals for personalized training plans")
    print("- These results are for educational purposes only")


def main(genome_file: str = "Genome.txt", output_json: Optional[str] = None):
    """
    Main function to run fitness analysis.
    
    Args:
        genome_file: Path to genome file
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze fitness
        results = analyze_fitness(genome_data)
        
        # Print report
        print_fitness_report(results)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'fitness': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Fitness & Athletic Performance Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, output_json=args.json)

