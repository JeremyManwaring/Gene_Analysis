"""Tests for probability calculation explanation."""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

import Probability_Calculation_Explanation as pce


def test_explain_probability_calculation():
    """Test probability calculation explanation."""
    result = pce.explain_probability_calculation()
    assert result == {
        'polygenic_score': 0.024,
        'aging_rate_probability': 0.976,
        'longevity_probability': 0.024,
        'cellular_health_probability': 0.024,
    }
