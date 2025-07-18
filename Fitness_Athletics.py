"""Utility wrapper exposing ``calculate_probability`` for tests."""


def calculate_probability(prs_score, baseline_prob=0.5):
    """Convert a PRS score to a probability percentage."""
    normalized_score = max(0, min(1, (prs_score + 1) / 2))
    probability = normalized_score * 100
    return round(probability, 1)
