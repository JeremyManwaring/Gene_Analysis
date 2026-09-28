"""Tests for fitness probability calculations."""

import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project_root))

from utils import calculate_probability


def test_calculate_probability_zero():
    """Test probability calculation with zero score."""
    assert calculate_probability(0) == 50.0


def test_calculate_probability_negative_one():
    """Test probability calculation with negative score."""
    assert calculate_probability(-1) == 0.0


def test_calculate_probability_one():
    """Test probability calculation with positive score."""
    assert calculate_probability(1) == 100.0
