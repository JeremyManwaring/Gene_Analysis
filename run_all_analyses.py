#!/usr/bin/env python3
"""
Genomic Analysis Runner
Runs disease and ethnicity analyses using shared functions.
"""

from disease_analysis import risk_summary, disease_statistics, polygenic_risk_scores
import subprocess
import sys


def run_ethnicity_analysis():
    """Run the Ethnicity analysis script (unchanged)."""
    try:
        result = subprocess.run([sys.executable, "Ethnicity.py"], capture_output=True, text=True, check=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running Ethnicity.py: {e}")
        print(f"Error output: {e.stderr}")


if __name__ == "__main__":
    print("GENOMIC SEQUENCE ANALYSIS")
    print("Analyzing your genetic data...")

    print(f"\n{'='*60}")
    print("RUNNING DISEASE ANALYSIS")
    print(f"{'='*60}")
    risk_summary()
    disease_statistics()
    polygenic_risk_scores()

    print(f"\n{'='*60}")
    print("RUNNING ETHNICITY ANALYSIS")
    print(f"{'='*60}")
    run_ethnicity_analysis()

    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print("Check the results above for your genetic insights!")
