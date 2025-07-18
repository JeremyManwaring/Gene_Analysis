import sys, os; sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import os
import pytest
from Longevity_Aging import LongevityAgingAnalyzer


SAMPLE_GENOME = {
    'rs2736100': 'GG',
    'rs7726159': 'AG',
    'rs1317082': 'AA',
    'rs1801133': 'AG',
    'rs1042522': 'AA',
}


def test_load_genome_data(tmp_path):
    sample_file = tmp_path / 'sample.txt'
    with open(sample_file, 'w') as f:
        for snp, gt in SAMPLE_GENOME.items():
            f.write(f"{snp} {gt}\n")
    analyzer = LongevityAgingAnalyzer()
    data = analyzer.load_genome_data(str(sample_file))
    assert data == SAMPLE_GENOME


def test_analyze_telomere_length():
    analyzer = LongevityAgingAnalyzer()
    result = analyzer.analyze_telomere_length(SAMPLE_GENOME)
    assert result['markers_found'] == 3
    assert result['total_markers'] == len(analyzer.telomere_markers)
    assert result['telomere_score'] == 0.21
    assert result['risk_level'] == 'High'


def test_calculate_polygenic_risk_score():
    analyzer = LongevityAgingAnalyzer()
    result = analyzer.calculate_polygenic_risk_score(SAMPLE_GENOME)
    assert result['markers_used'] == 5
    assert pytest.approx(result['polygenic_score'], rel=1e-6) == 0.305 / 26
    assert result['overall_aging_risk'] == 'High'
    assert pytest.approx(result['aging_rate_probability'], rel=1e-6) == 1 - (0.305 / 26)

