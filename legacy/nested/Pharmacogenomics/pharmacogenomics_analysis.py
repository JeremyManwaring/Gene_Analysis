if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Pharmacogenomics analysis - drug metabolism and response prediction."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data
from config.pharmacogenomics_snps import PHARMACOGENOMICS_SNPS, DRUG_RECOMMENDATIONS
from typing import Dict, Any, Optional, List


def analyze_pharmacogenomics(genome_data) -> Dict[str, Any]:
    """
    Analyze pharmacogenomics markers from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing pharmacogenomics analysis results
    """
    results = {
        'gene_variants': {},
        'drug_recommendations': {},
        'summary': {
            'total_markers': len(PHARMACOGENOMICS_SNPS),
            'markers_found': 0,
            'variants_detected': 0
        }
    }
    
    # Group by gene
    genes = {}
    for rsid, data in PHARMACOGENOMICS_SNPS.items():
        gene = data['gene']
        if gene not in genes:
            genes[gene] = []
        genes[gene].append((rsid, data))
    
    # Analyze each gene
    for gene, variants in genes.items():
        gene_results = {
            'gene': gene,
            'variants': [],
            'phenotype': 'Unknown',
            'drugs_affected': set()
        }
        
        for rsid, variant_data in variants:
            genotype = genome_data.get_genotype(rsid)
            
            variant_result = {
                'rsid': rsid,
                'variant': variant_data['variant'],
                'found': genotype is not None,
            }
            
            if genotype:
                results['summary']['markers_found'] += 1
                risk_allele = variant_data['risk_allele']
                has_variant = risk_allele in genotype
                
                variant_result['genotype'] = genotype
                variant_result['risk_allele'] = risk_allele
                variant_result['variant_present'] = has_variant
                variant_result['description'] = variant_data['description']
                
                if has_variant:
                    results['summary']['variants_detected'] += 1
                    # Determine if homozygous or heterozygous
                    allele_count = genotype.count(risk_allele)
                    variant_result['zygosity'] = 'Homozygous' if allele_count == 2 else 'Heterozygous'
                    variant_result['phenotype'] = variant_data.get('phenotype', 'Unknown')
                    
                    # Update gene phenotype (most significant variant)
                    if variant_data.get('phenotype'):
                        gene_results['phenotype'] = variant_data['phenotype']
                    
                    # Add affected drugs
                    for drug in variant_data.get('drugs', []):
                        gene_results['drugs_affected'].add(drug)
            else:
                variant_result['genotype'] = None
                variant_result['variant_present'] = False
            
            gene_results['variants'].append(variant_result)
        
        gene_results['drugs_affected'] = sorted(list(gene_results['drugs_affected']))
        results['gene_variants'][gene] = gene_results
    
    # Generate drug-specific recommendations
    for drug, drug_info in DRUG_RECOMMENDATIONS.items():
        affected_genes = drug_info['genes']
        recommendations = []
        
        for gene in affected_genes:
            if gene in results['gene_variants']:
                gene_data = results['gene_variants'][gene]
                if gene_data['phenotype'] != 'Unknown':
                    recommendations.append(f"{gene}: {gene_data['phenotype']}")
        
        if recommendations:
            results['drug_recommendations'][drug] = {
                'genes_analyzed': affected_genes,
                'findings': recommendations,
                'recommendation': drug_info['recommendation']
            }
    
    return results


def print_pharmacogenomics_report(results: Dict[str, Any]):
    """Print formatted pharmacogenomics report to console."""
    print("\n" + "="*60)
    print("PHARMACOGENOMICS ANALYSIS")
    print("="*60 + "\n")
    print("Drug Metabolism & Response Prediction")
    print("-"*60 + "\n")
    
    # Print gene variants
    for gene, gene_data in results['gene_variants'].items():
        print("="*60)
        print(f"GENE: {gene}")
        print("="*60)
        print(f"Phenotype: {gene_data['phenotype']}")
        if gene_data['drugs_affected']:
            print(f"Drugs Affected: {', '.join(gene_data['drugs_affected'])}")
        print()
        
        for variant in gene_data['variants']:
            print(f"  Variant: {variant['variant']}")
            if variant['found']:
                print(f"    Genotype: {variant['genotype']}")
                print(f"    Variant Present: {variant['variant_present']}")
                if variant['variant_present']:
                    print(f"    Zygosity: {variant.get('zygosity', 'N/A')}")
                    print(f"    Phenotype: {variant.get('phenotype', 'N/A')}")
                print(f"    Description: {variant['description']}")
            else:
                print(f"    Status: NOT FOUND IN DATA")
            print()
    
    # Print drug recommendations
    if results['drug_recommendations']:
        print("="*60)
        print("DRUG-SPECIFIC RECOMMENDATIONS")
        print("="*60 + "\n")
        
        for drug, rec_data in results['drug_recommendations'].items():
            print(f"Drug: {drug}")
            print(f"Genes Analyzed: {', '.join(rec_data['genes_analyzed'])}")
            print("Findings:")
            for finding in rec_data['findings']:
                print(f"  - {finding}")
            print(f"Recommendation: {rec_data['recommendation']}")
            print()
    
    # Summary
    summary = results['summary']
    print("="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total markers analyzed: {summary['total_markers']}")
    print(f"Markers found: {summary['markers_found']}")
    print(f"Variants detected: {summary['variants_detected']}")
    
    # Important disclaimer
    print("\n" + "="*60)
    print("IMPORTANT DISCLAIMER")
    print("="*60)
    print("- This analysis provides genetic information only")
    print("- Always consult healthcare providers before changing medications")
    print("- Genetic variants are one factor in drug response")
    print("- Other factors (age, health, other medications) also matter")
    print("- This is for educational purposes only, not medical advice")


def main(genome_file: str = "Genome.txt", output_json: Optional[str] = None):
    """
    Main function to run pharmacogenomics analysis.
    
    Args:
        genome_file: Path to genome file
        output_json: Optional path to save JSON results
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze pharmacogenomics
        results = analyze_pharmacogenomics(genome_data)
        
        # Print report
        print_pharmacogenomics_report(results)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'pharmacogenomics': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Pharmacogenomics Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--json", help="Export results to JSON file")
    
    args = parser.parse_args()
    main(genome_file=args.genome, output_json=args.json)

