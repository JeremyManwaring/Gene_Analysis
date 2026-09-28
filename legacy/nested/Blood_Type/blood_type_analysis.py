if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Blood type analysis (ABO and RHD systems)."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data
from config.blood_type_snps import ABO_SNPS, RHD_SNPS, BLOOD_TYPE_COMPATIBILITY
from typing import Dict, Any, Optional


def determine_abo_blood_type(genome_data) -> str:
    """
    Determine ABO blood type from genotype data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        ABO blood type (O, A, B, or AB)
    """
    # Check rs8176719 (O allele deletion)
    o_allele = genome_data.get_genotype('rs8176719')
    b_allele = genome_data.get_genotype('rs8176746')
    a_allele = genome_data.get_genotype('rs8176747')
    
    # Simplified logic - in reality, ABO typing is more complex
    # This is a simplified version for demonstration
    has_a = False
    has_b = False
    has_o = False
    
    if a_allele:
        # Check if A allele present (simplified)
        has_a = True
    
    if b_allele:
        # Check if B allele present (simplified)
        has_b = True
    
    if o_allele:
        # O allele deletion check (simplified)
        # In reality, need to check for deletion
        has_o = True
    
    # Determine blood type
    if has_a and has_b:
        return 'AB'
    elif has_a:
        return 'A'
    elif has_b:
        return 'B'
    else:
        return 'O'


def determine_rh_factor(genome_data) -> str:
    """
    Determine Rh factor from genotype data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Rh factor (+ or -)
    """
    rhd_genotype = genome_data.get_genotype('rs590787')
    
    if rhd_genotype:
        # Simplified: CC or CT = Rh+, TT = Rh-
        if 'T' in rhd_genotype and rhd_genotype.count('T') == 2:
            return '-'
        else:
            return '+'
    
    return 'Unknown'


def analyze_blood_type(genome_data) -> Dict[str, Any]:
    """
    Analyze blood type from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing blood type analysis results
    """
    results = {
        'abo_system': {},
        'rhd_system': {},
        'blood_type': 'Unknown',
        'compatibility': {}
    }
    
    # Analyze ABO system
    for rsid, data in ABO_SNPS.items():
        genotype = genome_data.get_genotype(rsid)
        results['abo_system'][rsid] = {
            'allele': data['allele'],
            'genotype': genotype,
            'found': genotype is not None
        }
    
    # Analyze RHD system
    for rsid, data in RHD_SNPS.items():
        genotype = genome_data.get_genotype(rsid)
        results['rhd_system'][rsid] = {
            'allele': data['allele'],
            'genotype': genotype,
            'found': genotype is not None
        }
    
    # Determine blood type
    abo_type = determine_abo_blood_type(genome_data)
    rh_type = determine_rh_factor(genome_data)
    
    if abo_type != 'Unknown' and rh_type != 'Unknown':
        results['blood_type'] = f"{abo_type}{rh_type}"
        
        # Get compatibility information
        if results['blood_type'] in BLOOD_TYPE_COMPATIBILITY:
            results['compatibility'] = BLOOD_TYPE_COMPATIBILITY[results['blood_type']]
    
    return results


def print_blood_type_report(results: Dict[str, Any]):
    """Print formatted blood type report to console."""
    print("\n" + "="*60)
    print("BLOOD TYPE ANALYSIS")
    print("="*60 + "\n")
    
    # Print blood type
    print("="*60)
    print("PREDICTED BLOOD TYPE")
    print("="*60)
    print(f"Blood Type: {results['blood_type']}")
    print()
    
    # Print ABO system
    print("="*60)
    print("ABO SYSTEM")
    print("="*60)
    for rsid, data in results['abo_system'].items():
        print(f"SNP: {rsid}")
        print(f"  Allele: {data['allele']}")
        if data['found']:
            print(f"  Genotype: {data['genotype']}")
        else:
            print(f"  Status: NOT FOUND IN DATA")
        print()
    
    # Print RHD system
    print("="*60)
    print("RHD SYSTEM (Rh Factor)")
    print("="*60)
    for rsid, data in results['rhd_system'].items():
        print(f"SNP: {rsid}")
        print(f"  Allele: {data['allele']}")
        if data['found']:
            print(f"  Genotype: {data['genotype']}")
        else:
            print(f"  Status: NOT FOUND IN DATA")
        print()
    
    # Print compatibility
    if results['compatibility']:
        print("="*60)
        print("BLOOD TYPE COMPATIBILITY")
        print("="*60)
        print(f"Can Receive From: {', '.join(results['compatibility'].get('can_receive', []))}")
        print(f"Can Donate To: {', '.join(results['compatibility'].get('can_donate_to', []))}")
        print()
    
    print("="*60)
    print("IMPORTANT NOTES")
    print("="*60)
    print("- This is a genetic prediction, not a lab test")
    print("- Always confirm blood type with medical laboratory testing")
    print("- Critical for blood transfusions and medical procedures")
    print("- This is for educational purposes only, not medical advice")


def main(genome_file: str = "Genome.txt", output_json: Optional[str] = None):
    """
    Main function to run blood type analysis.
    
    Args:
        genome_file: Path to genome file
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze blood type
        results = analyze_blood_type(genome_data)
        
        # Print report
        print_blood_type_report(results)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'blood_type': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Blood Type Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, output_json=args.json)

