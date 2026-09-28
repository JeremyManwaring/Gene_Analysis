"""Shared genotype calculation utilities."""

from typing import Dict, List, Tuple, Optional


def get_allele_count(genotype: str, risk_allele: str) -> int:
    """
    Count the number of risk alleles in a genotype string.
    
    Args:
        genotype: Genotype string (e.g., "AA", "AG", "GG")
        risk_allele: Risk allele to count (e.g., "A")
    
    Returns:
        Number of risk alleles (0, 1, or 2)
    """
    if not genotype or not risk_allele:
        return 0
    return genotype.count(risk_allele)


def genotype_effect(genotype: str, risk_allele: str) -> float:
    """
    Calculate effect score based on genotype and risk allele.
    
    Scoring:
    - Risk allele homozygous (e.g., AA when A is risk): 0.0
    - Heterozygous (e.g., AG when A is risk): 0.5
    - Protective homozygous (e.g., GG when A is risk): 1.0
    - Other/unknown: 0.3
    
    Args:
        genotype: Genotype string (e.g., "AA", "AG", "GG")
        risk_allele: Risk allele (e.g., "A")
    
    Returns:
        Effect score between 0.0 and 1.0
    """
    if not genotype or len(genotype) != 2:
        return 0.3
    
    # Risk allele homozygous
    if genotype == risk_allele * 2:
        return 0.0
    
    # Heterozygous
    if risk_allele in genotype:
        return 0.5
    
    # Determine protective allele (opposite of risk)
    # Common alleles: A, C, G, T
    protective_allele = _get_protective_allele(genotype, risk_allele)
    
    # Protective homozygous
    if protective_allele and genotype == protective_allele * 2:
        return 1.0
    
    # Other/unknown
    return 0.3


def _get_protective_allele(genotype: str, risk_allele: str) -> Optional[str]:
    """Determine protective allele from genotype and risk allele."""
    # Extract the non-risk allele
    alleles = set(genotype)
    if risk_allele in alleles:
        alleles.discard(risk_allele)
        if alleles:
            return list(alleles)[0]
    # If no risk allele found, return first allele as protective
    return genotype[0] if genotype else None


def calculate_prs(
    genotype_data: Dict[str, str],
    snp_weights: Dict[str, Tuple[str, float]],
    use_effect: bool = False
) -> Tuple[float, int]:
    """
    Calculate polygenic risk score from genotype data and SNP weights.
    
    Args:
        genotype_data: Dictionary mapping rsid to genotype
        snp_weights: Dictionary mapping rsid to (risk_allele, weight) tuple
        use_effect: If True, use genotype_effect() instead of allele count
    
    Returns:
        Tuple of (PRS score, number of SNPs found)
    """
    prs = 0.0
    snps_found = 0
    
    for rsid, (risk_allele, weight) in snp_weights.items():
        if rsid in genotype_data:
            genotype = genotype_data[rsid]
            if use_effect:
                count = genotype_effect(genotype, risk_allele)
            else:
                count = get_allele_count(genotype, risk_allele)
            prs += count * weight
            snps_found += 1
    
    return prs, snps_found


def normalize_score(score: float, total_markers: int, found_markers: Optional[int] = None) -> float:
    """
    Normalize a score by total markers.
    
    Args:
        score: Raw score to normalize
        total_markers: Total number of markers in the model
        found_markers: Number of markers actually found (if None, uses total_markers)
    
    Returns:
        Normalized score
    """
    if total_markers == 0:
        return 0.0
    
    if found_markers is None:
        found_markers = total_markers
    
    # Normalize by total markers (not found markers) for consistency
    return score / total_markers


def calculate_probability(prs_score: float, baseline_prob: float = 0.5) -> float:
    """
    Convert PRS score to probability percentage (0-100%).
    
    Args:
        prs_score: Polygenic risk score
        baseline_prob: Baseline probability (default: 0.5)
    
    Returns:
        Probability as percentage (0-100)
    """
    # Normalize PRS to probability scale
    normalized_score = max(0, min(1, (prs_score + 1) / 2))  # Scale to 0-1
    probability = normalized_score * 100
    return round(probability, 1)

