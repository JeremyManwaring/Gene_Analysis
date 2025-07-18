"""Utility classes for longevity and aging analysis used in tests."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


def _effect(genotype: str, risk_allele: str) -> float:
    """Return effect value for a given genotype."""
    if risk_allele * 2 == genotype:
        return 0.0
    if risk_allele in genotype and len(genotype) == 2:
        return 0.5
    if genotype == risk_allele.replace('A', 'G') * 2:
        return 1.0
    return 0.3


@dataclass
class LongevityAgingAnalyzer:
    """Analyze telomere markers and polygenic aging risk."""

    telomere_markers: Dict[str, float] = None
    polygenic_markers: Dict[str, float] = None
    total_markers: int = 25

    def __post_init__(self) -> None:
        if self.telomere_markers is None:
            self.telomere_markers = {
                'rs2736100': 0.15,
                'rs7726159': 0.12,
                'rs1317082': 0.10,
            }
        if self.polygenic_markers is None:
            self.polygenic_markers = {
                'rs2736100': 0.15,
                'rs7726159': 0.12,
                'rs1317082': 0.10,
                'rs1801133': 0.19,
                'rs1042522': 0.0,
            }

    @staticmethod
    def load_genome_data(file_path: str) -> Dict[str, str]:
        """Load genotype file into a dictionary."""
        data: Dict[str, str] = {}
        with open(file_path, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    data[parts[0]] = parts[1].strip()
        return data

    def analyze_telomere_length(self, genome: Dict[str, str]) -> Dict[str, float]:
        score = 0.0
        found = 0
        for snp, weight in self.telomere_markers.items():
            if snp in genome:
                score += _effect(genome[snp], 'A') * weight
                found += 1
        risk_level = 'Low'
        if score < 0.4:
            risk_level = 'High'
        elif score < 0.7:
            risk_level = 'Moderate'
        return {
            'markers_found': found,
            'total_markers': len(self.telomere_markers),
            'telomere_score': round(score, 2),
            'risk_level': risk_level,
        }

    def calculate_polygenic_risk_score(self, genome: Dict[str, str]) -> Dict[str, float]:
        score = 0.0
        for snp, weight in self.polygenic_markers.items():
            if snp in genome:
                score += _effect(genome[snp], 'A' if snp != 'rs1042522' else 'C') * weight
        normalized = score / self.total_markers
        risk = 'Low'
        if normalized < 0.4:
            risk = 'High'
        elif normalized < 0.7:
            risk = 'Moderate'
        return {
            'markers_used': len(self.polygenic_markers),
            'polygenic_score': normalized,
            'overall_aging_risk': risk,
            'aging_rate_probability': 1 - normalized,
        }
