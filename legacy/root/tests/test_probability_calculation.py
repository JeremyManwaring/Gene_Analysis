import sys, os; sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import Probability_Calculation_Explanation as pce


def test_explain_probability_calculation():
    result = pce.explain_probability_calculation()
    assert result == {
        'polygenic_score': 0.024,
        'aging_rate_probability': 0.976,
        'longevity_probability': 0.024,
        'cellular_health_probability': 0.024,
    }

