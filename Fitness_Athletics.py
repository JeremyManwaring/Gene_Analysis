import pandas as pd
import random

# Load genotype data
df = pd.read_csv('Genome.txt', sep='\t', comment='#', header=None)
df.columns = ['rsid', 'chromosome', 'position', 'genotype']

def get_allele_count(genotype, risk_allele):
    """Count the number of risk alleles in a genotype string."""
    return genotype.count(risk_allele)

def calculate_probability(prs_score, baseline_prob=0.5):
    """Convert PRS score to probability (0-100%)"""
    # Normalize PRS to probability scale
    normalized_score = max(0, min(1, (prs_score + 1) / 2))  # Scale to 0-1
    probability = normalized_score * 100
    return round(probability, 1)

# Fitness and Athletic Performance SNPs
fitness_snps = {
    'Muscle Fiber Type': {
        'rs1815739': {'risk': ['C'], 'description': 'ACTN3 - Fast twitch muscle fibers', 'weight': 0.30},
        'rs540874': {'risk': ['A'], 'description': 'ACE - Endurance performance', 'weight': 0.25},
        'rs4343': {'risk': ['G'], 'description': 'ACE - Power performance', 'weight': 0.20},
        'rs699': {'risk': ['A'], 'description': 'AGT - Blood pressure regulation', 'weight': 0.15},
        'rs5186': {'risk': ['C'], 'description': 'AGTR1 - Muscle hypertrophy', 'weight': 0.10},
    },
    'Endurance Capacity': {
        'rs1815739': {'risk': ['T'], 'description': 'ACTN3 - Endurance advantage', 'weight': 0.35},
        'rs540874': {'risk': ['A'], 'description': 'ACE - Endurance performance', 'weight': 0.30},
        'rs4343': {'risk': ['A'], 'description': 'ACE - Endurance variant', 'weight': 0.20},
        'rs2070744': {'risk': ['C'], 'description': 'NOS3 - Nitric oxide production', 'weight': 0.15},
    },
    'Power & Strength': {
        'rs1815739': {'risk': ['C'], 'description': 'ACTN3 - Power performance', 'weight': 0.40},
        'rs4343': {'risk': ['G'], 'description': 'ACE - Power variant', 'weight': 0.30},
        'rs5186': {'risk': ['C'], 'description': 'AGTR1 - Muscle strength', 'weight': 0.20},
        'rs699': {'risk': ['G'], 'description': 'AGT - Strength performance', 'weight': 0.10},
    },
    'Injury Risk': {
        'rs1800012': {'risk': ['A'], 'description': 'COL1A1 - ACL injury risk', 'weight': 0.35},
        'rs12722': {'risk': ['T'], 'description': 'COL5A1 - Tendon injury risk', 'weight': 0.30},
        'rs13946': {'risk': ['C'], 'description': 'COL12A1 - Ligament injury', 'weight': 0.20},
        'rs1800013': {'risk': ['G'], 'description': 'COL1A1 - Stress fracture risk', 'weight': 0.15},
        'rs2228570': {'risk': ['A'], 'description': 'VDR - Bone density and injury recovery', 'weight': 0.10},
    },
    'Recovery Rate': {
        'rs1800629': {'risk': ['A'], 'description': 'TNF - Inflammation response', 'weight': 0.30},
        'rs1801133': {'risk': ['T'], 'description': 'MTHFR - Recovery factors', 'weight': 0.25},
        'rs1801131': {'risk': ['C'], 'description': 'MTHFR - Recovery metabolism', 'weight': 0.20},
        'rs2070744': {'risk': ['T'], 'description': 'NOS3 - Recovery efficiency', 'weight': 0.15},
        'rs1799983': {'risk': ['T'], 'description': 'NOS3 - Recovery optimization', 'weight': 0.10},
    },
    'VO2 Max Potential': {
        'rs1815739': {'risk': ['T'], 'description': 'ACTN3 - Aerobic capacity', 'weight': 0.35},
        'rs540874': {'risk': ['A'], 'description': 'ACE - VO2 max', 'weight': 0.30},
        'rs4343': {'risk': ['A'], 'description': 'ACE - Endurance capacity', 'weight': 0.20},
        'rs2070744': {'risk': ['C'], 'description': 'NOS3 - Oxygen utilization', 'weight': 0.15},
    },
    'Muscle Hypertrophy': {
        'rs1815739': {'risk': ['C'], 'description': 'ACTN3 - Muscle growth', 'weight': 0.30},
        'rs5186': {'risk': ['C'], 'description': 'AGTR1 - Muscle hypertrophy', 'weight': 0.25},
        'rs4343': {'risk': ['G'], 'description': 'ACE - Muscle development', 'weight': 0.20},
        'rs699': {'risk': ['G'], 'description': 'AGT - Muscle building', 'weight': 0.15},
        'rs1800629': {'risk': ['G'], 'description': 'TNF - Muscle response', 'weight': 0.10},
    },
    'Lactate Threshold': {
        'rs1815739': {'risk': ['T'], 'description': 'ACTN3 - Lactate clearance', 'weight': 0.35},
        'rs540874': {'risk': ['A'], 'description': 'ACE - Lactate threshold', 'weight': 0.30},
        'rs4343': {'risk': ['A'], 'description': 'ACE - Endurance efficiency', 'weight': 0.20},
        'rs2070744': {'risk': ['C'], 'description': 'NOS3 - Metabolic efficiency', 'weight': 0.15},
    },
    'Caffeine Response': {
        'rs762551': {'risk': ['A'], 'description': 'CYP1A2 - Caffeine metabolism', 'weight': 0.40},
        'rs5751876': {'risk': ['T'], 'description': 'ADORA2A - Caffeine sensitivity', 'weight': 0.30},
        'rs5751879': {'risk': ['C'], 'description': 'ADORA2A - Performance response', 'weight': 0.20},
        'rs2470893': {'risk': ['T'], 'description': 'CYP1A2 - Caffeine effects', 'weight': 0.10},
    },
    'Heat Tolerance': {
        'rs1800629': {'risk': ['A'], 'description': 'TNF - Heat stress response', 'weight': 0.35},
        'rs1801133': {'risk': ['T'], 'description': 'MTHFR - Heat adaptation', 'weight': 0.25},
        'rs2070744': {'risk': ['C'], 'description': 'NOS3 - Thermoregulation', 'weight': 0.20},
        'rs1799983': {'risk': ['T'], 'description': 'NOS3 - Heat tolerance', 'weight': 0.20},
    }
}

print("\n==============================")
print("FITNESS & ATHLETIC PERFORMANCE ANALYSIS")
print("==============================\n")

# Track results for summary
all_results = {}
trait_probabilities = {}

for trait, snps in fitness_snps.items():
    print(f"{'='*60}")
    print(f"TRAIT: {trait}")
    print(f"{'='*60}")
    
    prs = 0.0
    snps_found = 0
    snp_details = []
    
    for rsid, data in snps.items():
        match = df[df['rsid'] == rsid]
        if not match.empty:
            genotype = match.iloc[0]['genotype']
            count = get_allele_count(genotype, data['risk'][0])
            prs += count * data['weight']
            snps_found += 1
            snp_details.append((rsid, genotype, data['risk'][0], data['weight'], count, data['description']))
            all_results[rsid] = (genotype, data['description'], trait)
        else:
            snp_details.append((rsid, 'Not found', data['risk'][0], data['weight'], 0, data['description']))
    
    # Print SNP table
    print(f"{'SNP':<12}{'Genotype':<12}{'Risk':<8}{'Weight':<10}{'Count':<8}{'Description':<30}")
    print(f"{'-'*85}")
    for rsid, genotype, risk_allele, weight, count, description in snp_details:
        print(f"{rsid:<12}{genotype:<12}{risk_allele:<8}{weight:<10}{count:<8}{description:<30}")
    print(f"{'-'*85}")
    
    # Calculate probability
    probability = calculate_probability(prs)
    trait_probabilities[trait] = probability
    
    print(f"PRS for {trait}: {prs:.2f} (based on {snps_found}/{len(snps)} SNPs found)")
    print(f"Probability Score: {probability}%")
    print()

# Summary and Recommendations
print(f"{'='*60}")
print("SUMMARY & PROBABILITY SCORES")
print(f"{'='*60}")

print(f"\n{'Trait':<25}{'Probability':<15}{'Interpretation'}")
print(f"{'-'*60}")

for trait, probability in trait_probabilities.items():
    if probability >= 70:
        interpretation = "HIGH"
    elif probability >= 50:
        interpretation = "MODERATE"
    else:
        interpretation = "LOW"
    
    print(f"{trait:<25}{probability:<15}%{interpretation}")

print(f"\n{'='*60}")
print("DETAILED INTERPRETATION")
print(f"{'='*60}")

for trait, probability in trait_probabilities.items():
    print(f"\n{trait}:")
    if trait == 'Muscle Fiber Type':
        if probability >= 70:
            print(f"  {probability}% probability of FAST-TWITCH muscle dominance")
            print("  → Better suited for power sports (sprinting, weightlifting)")
        elif probability >= 50:
            print(f"  {probability}% probability of MIXED muscle fiber type")
            print("  → Balanced for both power and endurance")
        else:
            print(f"  {probability}% probability of SLOW-TWITCH muscle dominance")
            print("  → Better suited for endurance sports (distance running, cycling)")
    
    elif trait == 'Endurance Capacity':
        if probability >= 70:
            print(f"  {probability}% probability of HIGH endurance capacity")
            print("  → Excellent for long-distance events")
        elif probability >= 50:
            print(f"  {probability}% probability of MODERATE endurance capacity")
            print("  → Good baseline for endurance training")
        else:
            print(f"  {probability}% probability of LOWER endurance capacity")
            print("  → May need more training for endurance events")
    
    elif trait == 'Power & Strength':
        if probability >= 70:
            print(f"  {probability}% probability of HIGH power potential")
            print("  → Excellent for explosive movements")
        elif probability >= 50:
            print(f"  {probability}% probability of MODERATE power potential")
            print("  → Good baseline for strength training")
        else:
            print(f"  {probability}% probability of LOWER power potential")
            print("  → May need more focus on strength training")
    
    elif trait == 'Injury Risk':
        if probability >= 70:
            print(f"  {probability}% probability of HIGHER injury risk")
            print("  → Focus on injury prevention and proper form")
        elif probability >= 50:
            print(f"  {probability}% probability of MODERATE injury risk")
            print("  → Standard injury prevention recommended")
        else:
            print(f"  {probability}% probability of LOWER injury risk")
            print("  → Good genetic resilience, but still practice safety")
    
    elif trait == 'Recovery Rate':
        if probability >= 70:
            print(f"  {probability}% probability of FAST recovery")
            print("  → Can handle higher training volumes")
        elif probability >= 50:
            print(f"  {probability}% probability of MODERATE recovery")
            print("  → Standard recovery protocols recommended")
        else:
            print(f"  {probability}% probability of SLOWER recovery")
            print("  → May need more recovery time between sessions")
    
    else:
        if probability >= 70:
            print(f"  {probability}% probability of HIGH {trait.lower()}")
        elif probability >= 50:
            print(f"  {probability}% probability of MODERATE {trait.lower()}")
        else:
            print(f"  {probability}% probability of LOWER {trait.lower()}")

print(f"\n{'='*60}")
print("TRAINING RECOMMENDATIONS")
print(f"{'='*60}")
print("- These probabilities indicate genetic predispositions, not limitations")
print("- Training can significantly improve performance regardless of genetics")
print("- Focus on your strengths while developing areas of opportunity")
print("- Consult with fitness professionals for personalized training plans")
print("- These results are for educational purposes only") 