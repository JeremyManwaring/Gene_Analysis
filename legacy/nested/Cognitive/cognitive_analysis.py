if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")

"""Comprehensive cognitive and intelligence analysis."""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from utils import load_genome_data, genotype_effect
from config.cognitive_snps import COGNITIVE_SNPS, COGNITIVE_DOMAIN_WEIGHTS
from typing import Dict, Any, Optional, List, Tuple


def analyze_cognitive_traits(genome_data) -> Dict[str, Any]:
    """
    Analyze cognitive traits from genome data.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing cognitive analysis results
    """
    results = {
        'cognitive_domains': {},
        'domain_scores': {},
        'overall_cognitive_score': 0.0,
        'summary': {
            'total_markers': sum(len(snps) for snps in COGNITIVE_SNPS.values()),
            'markers_found': 0,
            'variants_detected': 0
        }
    }
    
    total_weighted_score = 0.0
    total_weight = 0.0
    
    # Analyze each cognitive domain
    for domain, snps in COGNITIVE_SNPS.items():
        domain_results = {
            'domain': domain,
            'markers': [],
            'domain_score': 0.0,
            'markers_found': 0
        }
        
        domain_score = 0.0
        domain_weight = COGNITIVE_DOMAIN_WEIGHTS.get(domain, 0.1)
        
        for rsid, snp_data in snps.items():
            genotype = genome_data.get_genotype(rsid)
            
            marker_result = {
                'rsid': rsid,
                'gene': snp_data['gene'],
                'variant': snp_data['variant'],
                'found': genotype is not None,
            }
            
            if genotype:
                results['summary']['markers_found'] += 1
                domain_results['markers_found'] += 1
                
                risk_allele = snp_data['risk_allele']
                has_variant = risk_allele in genotype
                
                marker_result['genotype'] = genotype
                marker_result['risk_allele'] = risk_allele
                marker_result['variant_present'] = has_variant
                marker_result['description'] = snp_data['description']
                marker_result['cognitive_impact'] = snp_data.get('cognitive_impact', '')
                
                # Get interpretation based on genotype
                if 'interpretation' in snp_data:
                    marker_result['interpretation'] = snp_data['interpretation'].get(
                        genotype, 'Genotype interpretation not available'
                    )
                
                # Calculate effect score (inverse for cognitive decline risk)
                if domain == 'Cognitive Decline Risk':
                    # For decline risk, lower is better
                    effect = genotype_effect(genotype, risk_allele)
                    # Invert: 1.0 becomes 0.0 (good), 0.0 becomes 1.0 (bad)
                    effect = 1.0 - effect
                else:
                    # For other domains, higher is better
                    effect = genotype_effect(genotype, risk_allele)
                
                marker_result['effect_score'] = round(effect, 3)
                domain_score += effect
                
                if has_variant:
                    results['summary']['variants_detected'] += 1
            else:
                marker_result['genotype'] = None
                marker_result['variant_present'] = False
                marker_result['effect_score'] = 0.0
            
            domain_results['markers'].append(marker_result)
        
        # Normalize domain score
        if domain_results['markers_found'] > 0:
            domain_results['domain_score'] = round(
                domain_score / domain_results['markers_found'], 3
            )
            total_weighted_score += domain_results['domain_score'] * domain_weight
            total_weight += domain_weight
        
        results['cognitive_domains'][domain] = domain_results
        results['domain_scores'][domain] = domain_results['domain_score']
    
    # Calculate overall cognitive score
    if total_weight > 0:
        results['overall_cognitive_score'] = round(total_weighted_score / total_weight, 3)
    
    # Determine cognitive profile
    if results['overall_cognitive_score'] >= 0.7:
        results['cognitive_profile'] = 'High'
    elif results['overall_cognitive_score'] >= 0.5:
        results['cognitive_profile'] = 'Moderate'
    else:
        results['cognitive_profile'] = 'Standard'
    
    # Assess cognitive decline risk separately
    if 'Cognitive Decline Risk' in results['domain_scores']:
        decline_score = results['domain_scores']['Cognitive Decline Risk']
        if decline_score < 0.4:
            results['decline_risk'] = 'High'
        elif decline_score < 0.7:
            results['decline_risk'] = 'Moderate'
        else:
            results['decline_risk'] = 'Low'
    else:
        results['decline_risk'] = 'Unknown'
    
    return results


def analyze_learning_ability_detailed(genome_data) -> Dict[str, Any]:
    """
    Detailed analysis of learning ability markers with deep dive into combinations.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary containing detailed learning ability analysis
    """
    learning_snps = COGNITIVE_SNPS.get('Learning Ability', {})
    
    results = {
        'neuroplasticity': {},
        'reward_learning': {},
        'memory_consolidation': {},
        'attention_focus': {},
        'combined_analysis': {},
        'learning_style': {},
        'nutrition_impact': {}
    }
    
    # Analyze by category
    categories = {
        'Neuroplasticity': [],
        'Reward Learning': [],
        'Memory Consolidation': [],
        'Attention/Focus': []
    }
    
    for rsid, snp_data in learning_snps.items():
        category = snp_data.get('category', 'Unknown')
        genotype = genome_data.get_genotype(rsid)
        
        marker_result = {
            'rsid': rsid,
            'gene': snp_data['gene'],
            'variant': snp_data['variant'],
            'found': genotype is not None,
            'category': category
        }
        
        if genotype:
            risk_allele = snp_data['risk_allele']
            effect = genotype_effect(genotype, risk_allele)
            
            marker_result['genotype'] = genotype
            marker_result['effect_score'] = round(effect, 3)
            marker_result['description'] = snp_data['description']
            marker_result['cognitive_impact'] = snp_data.get('cognitive_impact', '')
            
            if 'interpretation' in snp_data:
                marker_result['interpretation'] = snp_data['interpretation'].get(
                    genotype, 'Genotype interpretation not available'
                )
            
            # Add to category
            if category in categories:
                categories[category].append(marker_result)
        
        # Store in appropriate result section
        if category == 'Neuroplasticity':
            results['neuroplasticity'][rsid] = marker_result
        elif category == 'Reward Learning':
            results['reward_learning'][rsid] = marker_result
        elif category == 'Memory Consolidation':
            results['memory_consolidation'][rsid] = marker_result
        elif category == 'Attention/Focus':
            results['attention_focus'][rsid] = marker_result
    
    # Calculate category scores
    for category, markers in categories.items():
        if markers:
            found_markers = [m for m in markers if m['found']]
            if found_markers:
                avg_score = sum(m['effect_score'] for m in found_markers) / len(found_markers)
                results['combined_analysis'][category.lower().replace('/', '_')] = {
                    'score': round(avg_score, 3),
                    'markers_found': len(found_markers),
                    'total_markers': len(markers)
                }
    
    # Analyze BDNF and DRD2 combination
    bdnf_genotype = genome_data.get_genotype('rs6265')
    drd2_genotype = genome_data.get_genotype('rs1800497')
    
    if bdnf_genotype and drd2_genotype:
        bdnf_effect = genotype_effect(bdnf_genotype, 'C')  # Met is risk
        drd2_effect = genotype_effect(drd2_genotype, 'A')  # A1 is risk
        
        # Combined effect interpretation
        if bdnf_effect >= 0.7 and drd2_effect >= 0.7:
            combination = 'Optimal Learning Profile'
            description = 'High neuroplasticity + High reward learning = Excellent learning potential'
        elif bdnf_effect >= 0.7 or drd2_effect >= 0.7:
            combination = 'Strong Learning Profile'
            description = 'One strong learning factor present'
        else:
            combination = 'Standard Learning Profile'
            description = 'Standard learning genetics'
        
        results['combined_analysis']['bdnf_drd2_combination'] = {
            'combination': combination,
            'description': description,
            'bdnf_score': round(bdnf_effect, 3),
            'drd2_score': round(drd2_effect, 3),
            'combined_score': round((bdnf_effect + drd2_effect) / 2, 3)
        }
    
    # Predict learning style
    results['learning_style'] = predict_learning_style(genome_data)
    
    # Analyze nutrition impact
    results['nutrition_impact'] = analyze_nutrition_impact_on_learning(genome_data)
    
    return results


def predict_learning_style(genome_data) -> Dict[str, Any]:
    """
    Predict learning style based on genetic markers.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary with predicted learning style and reasoning
    """
    scores = {
        'visual': 0.0,
        'auditory': 0.0,
        'kinesthetic': 0.0,
        'mixed': 0.0
    }
    
    reasoning = []
    
    # BDNF analysis (affects neuroplasticity - visual learning)
    bdnf_genotype = genome_data.get_genotype('rs6265')
    if bdnf_genotype:
        bdnf_effect = genotype_effect(bdnf_genotype, 'C')
        if bdnf_effect >= 0.7:
            scores['visual'] += 0.4
            reasoning.append('High BDNF (neuroplasticity) supports visual learning')
        scores['mixed'] += bdnf_effect * 0.3
    
    # DRD2 analysis (affects reward learning - auditory/feedback)
    drd2_genotype = genome_data.get_genotype('rs1800497')
    if drd2_genotype:
        drd2_effect = genotype_effect(drd2_genotype, 'A')
        if drd2_effect >= 0.7:
            scores['auditory'] += 0.4
            reasoning.append('High DRD2 (reward learning) supports auditory/feedback learning')
        scores['mixed'] += drd2_effect * 0.3
    
    # COMT analysis (affects working memory - kinesthetic)
    comt_genotype = genome_data.get_genotype('rs4680')
    if comt_genotype:
        comt_effect = genotype_effect(comt_genotype, 'G')
        if comt_effect <= 0.3:  # Val/Val (GG) - better under stress
            scores['kinesthetic'] += 0.3
            reasoning.append('COMT Val/Val supports kinesthetic/hands-on learning')
        scores['mixed'] += (1 - comt_effect) * 0.2
    
    # Determine primary learning style
    max_score = max(scores.values())
    primary_style = [style for style, score in scores.items() if score == max_score][0]
    
    # If scores are close, prefer mixed
    if max_score < 0.5 or abs(max(scores.values()) - min(scores.values())) < 0.2:
        primary_style = 'mixed'
    
    return {
        'primary_style': primary_style.capitalize(),
        'scores': {k: round(v, 3) for k, v in scores.items()},
        'reasoning': reasoning,
        'recommendations': get_learning_style_recommendations(primary_style)
    }


def get_learning_style_recommendations(style: str) -> List[str]:
    """Get learning recommendations based on predicted style."""
    recommendations = {
        'visual': [
            'Use visual aids: diagrams, charts, mind maps',
            'Color-code information',
            'Watch videos and demonstrations',
            'Create visual summaries',
            'Use flashcards with images'
        ],
        'auditory': [
            'Listen to lectures and podcasts',
            'Read aloud or use text-to-speech',
            'Participate in discussions',
            'Use mnemonic devices and rhymes',
            'Record and replay information'
        ],
        'kinesthetic': [
            'Hands-on activities and experiments',
            'Take notes while learning',
            'Use movement while studying',
            'Build models or prototypes',
            'Practice through doing'
        ],
        'mixed': [
            'Combine multiple learning methods',
            'Use visual, auditory, and hands-on approaches',
            'Adapt to different learning contexts',
            'Experiment with different techniques',
            'Use varied study strategies'
        ]
    }
    return recommendations.get(style.lower(), recommendations['mixed'])


def analyze_nutrition_impact_on_learning(genome_data) -> Dict[str, Any]:
    """
    Analyze how nutrition markers affect learning ability.
    
    Args:
        genome_data: GenomeData object with genotype information
    
    Returns:
        Dictionary with nutrition impact on learning
    """
    impact = {
        'factors': [],
        'recommendations': [],
        'combined_impact': 'Neutral'
    }
    
    # Check MTHFR (folate metabolism)
    mthfr_677 = genome_data.get_genotype('rs1801133')
    mthfr_1298 = genome_data.get_genotype('rs1801131')
    
    folate_impact = False
    if mthfr_677 or mthfr_1298:
        if mthfr_677 and ('T' in mthfr_677 or mthfr_677 == 'TT'):
            folate_impact = True
            impact['factors'].append({
                'nutrient': 'Folate',
                'variant': 'MTHFR C677T',
                'impact': 'Reduced folate conversion may affect BDNF and learning',
                'recommendation': 'Use methylfolate instead of folic acid'
            })
            impact['recommendations'].append('Consider methylfolate supplementation for optimal BDNF function')
        
        if mthfr_1298 and ('C' in mthfr_1298 or mthfr_1298 == 'CC'):
            folate_impact = True
            impact['factors'].append({
                'nutrient': 'Folate',
                'variant': 'MTHFR A1298C',
                'impact': 'Mildly reduced folate metabolism',
                'recommendation': 'Monitor folate levels'
            })
    
    # Check B12
    b12_genotype = genome_data.get_genotype('rs602662')
    if b12_genotype and 'G' in b12_genotype:
        impact['factors'].append({
            'nutrient': 'Vitamin B12',
            'variant': 'FUT2',
            'impact': 'Lower B12 levels may affect cognitive function',
            'recommendation': 'Monitor B12 levels, may need supplementation'
        })
        impact['recommendations'].append('Ensure adequate B12 for optimal learning and memory')
    
    # Check Vitamin D
    vitd_genotype = genome_data.get_genotype('rs2282679')
    if vitd_genotype and 'C' in vitd_genotype:
        impact['factors'].append({
            'nutrient': 'Vitamin D',
            'variant': 'GC',
            'impact': 'Lower vitamin D may affect BDNF and cognitive function',
            'recommendation': 'Consider vitamin D supplementation (2000-4000 IU)'
        })
        impact['recommendations'].append('Maintain optimal vitamin D levels for brain health and learning')
    
    # Determine combined impact
    if len(impact['factors']) >= 2:
        impact['combined_impact'] = 'Significant - Multiple nutrition factors may affect learning'
    elif len(impact['factors']) == 1:
        impact['combined_impact'] = 'Moderate - One nutrition factor may affect learning'
    else:
        impact['combined_impact'] = 'Minimal - Standard nutrition needs'
    
    return impact


def get_learning_recommendations(results: Dict[str, Any]) -> list:
    """Generate learning recommendations based on cognitive profile."""
    recommendations = []
    
    # Working memory recommendations
    if 'Working Memory' in results['domain_scores']:
        wm_score = results['domain_scores']['Working Memory']
        if wm_score >= 0.7:
            recommendations.append("Strong working memory - can handle complex multi-step tasks")
        elif wm_score < 0.5:
            recommendations.append("May benefit from breaking tasks into smaller steps")
    
    # Memory recommendations
    if 'Long-term Memory' in results['domain_scores']:
        mem_score = results['domain_scores']['Long-term Memory']
        if mem_score >= 0.7:
            recommendations.append("Strong memory - effective use of spaced repetition recommended")
        elif mem_score < 0.5:
            recommendations.append("May benefit from increased repetition and review")
    
    # Learning ability recommendations
    if 'Learning Ability' in results['domain_scores']:
        learn_score = results['domain_scores']['Learning Ability']
        if learn_score >= 0.7:
            recommendations.append("High learning ability - can adapt to new concepts quickly")
        elif learn_score < 0.5:
            recommendations.append("May benefit from slower-paced, structured learning")
    
    # Cognitive decline risk recommendations
    if results.get('decline_risk') == 'High':
        recommendations.append("Higher cognitive decline risk - maintain brain health activities")
        recommendations.append("Consider: regular exercise, mental stimulation, social engagement")
        recommendations.append("Monitor nutrition (folate, B12, omega-3) for brain health")
    
    return recommendations


def print_learning_ability_report(detailed_results: Dict[str, Any]):
    """Print detailed learning ability analysis report."""
    print("\n" + "="*60)
    print("DETAILED LEARNING ABILITY ANALYSIS")
    print("="*60 + "\n")
    
    # BDNF Analysis
    if 'rs6265' in detailed_results.get('neuroplasticity', {}):
        bdnf = detailed_results['neuroplasticity']['rs6265']
        print("="*60)
        print("BDNF (BRAIN-DERIVED NEUROTROPHIC FACTOR) ANALYSIS")
        print("="*60)
        if bdnf['found']:
            print(f"Genotype: {bdnf['genotype']}")
            print(f"Effect Score: {bdnf['effect_score']:.3f}")
            print(f"Interpretation: {bdnf.get('interpretation', 'N/A')}")
            print(f"\nWhat This Means:")
            print(f"  • BDNF is crucial for neuroplasticity - your brain's ability to adapt")
            print(f"  • {bdnf['genotype']} genotype indicates: {bdnf.get('interpretation', 'N/A')}")
            print(f"  • Higher BDNF = Better memory formation and learning capacity")
            print(f"  • Responds well to: Exercise, learning new skills, cognitive training")
        else:
            print("BDNF marker not found in data")
        print()
    
    # DRD2 Analysis
    if 'rs1800497' in detailed_results.get('reward_learning', {}):
        drd2 = detailed_results['reward_learning']['rs1800497']
        print("="*60)
        print("DRD2 (DOPAMINE RECEPTOR D2) ANALYSIS")
        print("="*60)
        if drd2['found']:
            print(f"Genotype: {drd2['genotype']}")
            print(f"Effect Score: {drd2['effect_score']:.3f}")
            print(f"Interpretation: {drd2.get('interpretation', 'N/A')}")
            print(f"\nWhat This Means:")
            print(f"  • DRD2 affects reward processing and motivation")
            print(f"  • {drd2['genotype']} genotype indicates: {drd2.get('interpretation', 'N/A')}")
            print(f"  • Higher DRD2 = Better learning from feedback and rewards")
            print(f"  • Responds well to: Reward-based learning, positive feedback, achievement systems")
        else:
            print("DRD2 marker not found in data")
        print()
    
    # Combined BDNF + DRD2 Analysis
    if 'bdnf_drd2_combination' in detailed_results.get('combined_analysis', {}):
        combo = detailed_results['combined_analysis']['bdnf_drd2_combination']
        print("="*60)
        print("COMBINED BDNF + DRD2 EFFECT")
        print("="*60)
        print(f"Combination Profile: {combo['combination']}")
        print(f"Description: {combo['description']}")
        print(f"BDNF Score: {combo['bdnf_score']:.3f}")
        print(f"DRD2 Score: {combo['drd2_score']:.3f}")
        print(f"Combined Score: {combo['combined_score']:.3f}")
        print()
        
        if combo['combined_score'] >= 0.7:
            print("🎯 OPTIMAL LEARNING PROFILE DETECTED")
            print("  • High neuroplasticity + High reward learning")
            print("  • Excellent potential for rapid skill acquisition")
            print("  • Responds exceptionally well to structured learning")
            print("  • Benefits from both visual and feedback-based methods")
        print()
    
    # Learning Style Prediction
    if 'learning_style' in detailed_results:
        style = detailed_results['learning_style']
        print("="*60)
        print("PREDICTED LEARNING STYLE")
        print("="*60)
        print(f"Primary Style: {style['primary_style']}")
        print(f"\nStyle Scores:")
        for s, score in style['scores'].items():
            print(f"  {s.capitalize()}: {score:.3f}")
        print(f"\nReasoning:")
        for reason in style.get('reasoning', []):
            print(f"  • {reason}")
        print(f"\nRecommended Learning Strategies:")
        for rec in style.get('recommendations', []):
            print(f"  • {rec}")
        print()
    
    # Category Analysis
    print("="*60)
    print("LEARNING ABILITY BY CATEGORY")
    print("="*60)
    for category, data in detailed_results.get('combined_analysis', {}).items():
        if isinstance(data, dict) and 'score' in data:
            print(f"{category.replace('_', ' ').title()}:")
            print(f"  Score: {data['score']:.3f}")
            print(f"  Markers Found: {data['markers_found']}/{data['total_markers']}")
            print()
    
    # Nutrition Impact
    if 'nutrition_impact' in detailed_results:
        nutrition = detailed_results['nutrition_impact']
        print("="*60)
        print("NUTRITION IMPACT ON LEARNING")
        print("="*60)
        print(f"Combined Impact: {nutrition['combined_impact']}")
        if nutrition['factors']:
            print(f"\nNutrition Factors Affecting Learning:")
            for factor in nutrition['factors']:
                print(f"  • {factor['nutrient']} ({factor['variant']})")
                print(f"    Impact: {factor['impact']}")
                print(f"    Recommendation: {factor['recommendation']}")
                print()
        if nutrition['recommendations']:
            print("Nutrition Recommendations for Optimal Learning:")
            for rec in nutrition['recommendations']:
                print(f"  • {rec}")
        print()
    
    # Optimal Learning Strategies
    print("="*60)
    print("OPTIMAL LEARNING STRATEGIES FOR YOUR GENETIC PROFILE")
    print("="*60)
    
    # Get BDNF and DRD2 scores
    bdnf_score = 0.0
    drd2_score = 0.0
    
    if 'rs6265' in detailed_results.get('neuroplasticity', {}):
        bdnf_score = detailed_results['neuroplasticity']['rs6265'].get('effect_score', 0.0)
    if 'rs1800497' in detailed_results.get('reward_learning', {}):
        drd2_score = detailed_results['reward_learning']['rs1800497'].get('effect_score', 0.0)
    
    if bdnf_score >= 0.7:
        print("\n✓ HIGH NEUROPLASTICITY PROFILE")
        print("  • Your brain adapts quickly to new information")
        print("  • Strategies:")
        print("    - Learn multiple skills simultaneously")
        print("    - Use spaced repetition (your brain consolidates well)")
        print("    - Engage in regular exercise (boosts BDNF further)")
        print("    - Challenge yourself with new, complex material")
    
    if drd2_score >= 0.7:
        print("\n✓ HIGH REWARD LEARNING PROFILE")
        print("  • You learn best from feedback and rewards")
        print("  • Strategies:")
        print("    - Use gamification and achievement systems")
        print("    - Set up reward milestones")
        print("    - Seek immediate feedback on performance")
        print("    - Track progress visually")
    
    if bdnf_score >= 0.7 and drd2_score >= 0.7:
        print("\n🎯 COMBINED ADVANTAGE")
        print("  • You have BOTH high neuroplasticity AND reward learning")
        print("  • This is an optimal combination for learning")
        print("  • You can:")
        print("    - Learn quickly AND stay motivated")
        print("    - Benefit from both structured and exploratory learning")
        print("    - Excel in environments with clear feedback and rewards")
    
    print("\n" + "="*60)
    print("ENVIRONMENTAL FACTORS FOR OPTIMAL LEARNING")
    print("="*60)
    print("• Sleep: Critical for memory consolidation (aim for 7-9 hours)")
    print("• Exercise: Regular physical activity boosts BDNF")
    print("• Stress Management: Your COMT profile helps, but still important")
    print("• Nutrition: Follow recommendations above for folate, B12, Vitamin D")
    print("• Social Learning: Interaction with others enhances BDNF")
    print("• Novelty: Learning new things regularly maintains neuroplasticity")
    print()


def print_cognitive_report(results: Dict[str, Any]):
    """Print formatted cognitive analysis report to console."""
    print("\n" + "="*60)
    print("COGNITIVE & INTELLIGENCE ANALYSIS")
    print("="*60 + "\n")
    
    # Run detailed learning ability analysis
    from utils import load_genome_data
    # Get genome data from results if available, otherwise we'll need to reload
    # For now, we'll add the detailed analysis at the end
    
    # Overall cognitive profile
    print("="*60)
    print("OVERALL COGNITIVE PROFILE")
    print("="*60)
    print(f"Overall Cognitive Score: {results['overall_cognitive_score']:.3f}")
    print(f"Profile: {results['cognitive_profile']}")
    print(f"Cognitive Decline Risk: {results.get('decline_risk', 'Unknown')}")
    print()
    
    # Domain scores summary
    print("="*60)
    print("COGNITIVE DOMAIN SCORES")
    print("="*60)
    print(f"{'Domain':<25} {'Score':<10} {'Interpretation'}")
    print("-"*60)
    
    for domain, score in results['domain_scores'].items():
        if score >= 0.7:
            interpretation = "High"
        elif score >= 0.5:
            interpretation = "Moderate"
        else:
            interpretation = "Standard"
        print(f"{domain:<25} {score:<10.3f} {interpretation}")
    print()
    
    # Detailed domain analysis
    for domain, domain_data in results['cognitive_domains'].items():
        print("="*60)
        print(f"DOMAIN: {domain}")
        print("="*60)
        print(f"Domain Score: {domain_data['domain_score']:.3f}")
        print(f"Markers Found: {domain_data['markers_found']}/{len(domain_data['markers'])}")
        print()
        
        for marker in domain_data['markers']:
            print(f"Gene: {marker['gene']} ({marker['rsid']})")
            print(f"  Variant: {marker['variant']}")
            if marker['found']:
                print(f"  Genotype: {marker['genotype']}")
                print(f"  Effect Score: {marker['effect_score']:.3f}")
                if marker.get('interpretation'):
                    print(f"  Interpretation: {marker['interpretation']}")
                print(f"  Description: {marker['description']}")
                print(f"  Cognitive Impact: {marker.get('cognitive_impact', 'N/A')}")
            else:
                print(f"  Status: NOT FOUND IN DATA")
            print()
    
    # Learning recommendations
    recommendations = get_learning_recommendations(results)
    if recommendations:
        print("="*60)
        print("LEARNING RECOMMENDATIONS")
        print("="*60)
        for rec in recommendations:
            print(f"- {rec}")
        print()
    
    # Summary
    summary = results['summary']
    print("="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Total markers analyzed: {summary['total_markers']}")
    print(f"Markers found: {summary['markers_found']}")
    print(f"Variants detected: {summary['variants_detected']}")
    
    # Important disclaimers
    print("\n" + "="*60)
    print("IMPORTANT DISCLAIMERS")
    print("="*60)
    print("- Intelligence is highly polygenic (involves thousands of genetic variants)")
    print("- Environment, education, and experience play major roles in cognitive ability")
    print("- These genetic markers show predispositions, not fixed abilities")
    print("- Cognitive performance can be improved through training and practice")
    print("- These results are for educational purposes only, not definitive IQ scores")
    print("- Consult educational professionals for learning strategies")


def main(genome_file: str = "Genome.txt", output_json: Optional[str] = None, detailed_learning: bool = True):
    """
    Main function to run cognitive analysis.
    
    Args:
        genome_file: Path to genome file
        output_json: Optional path to save JSON results
        detailed_learning: Include detailed learning ability analysis
    """
    try:
        # Load genome data
        genome_data = load_genome_data(genome_file)
        
        # Analyze cognitive traits
        results = analyze_cognitive_traits(genome_data)
        
        # Print main cognitive report
        print_cognitive_report(results)
        
        # Detailed learning ability analysis
        if detailed_learning:
            learning_detailed = analyze_learning_ability_detailed(genome_data)
            results['learning_ability_detailed'] = learning_detailed
            print_learning_ability_report(learning_detailed)
        
        # Export to JSON if requested
        if output_json:
            from utils import export_to_json
            export_to_json({'cognitive': results}, filename=output_json)
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
    
    parser = argparse.ArgumentParser(description="Cognitive & Intelligence Analysis")
    parser.add_argument("--genome", default="Genome.txt", help="Path to genome file")
    parser.add_argument("--json", help="Export results to JSON file")
    parser.add_argument("--no-detailed-learning", action="store_true", help="Skip detailed learning ability analysis")
    
    args = parser.parse_args()
    main(genome_file=args.genome, output_json=args.json, detailed_learning=not args.no_detailed_learning)

