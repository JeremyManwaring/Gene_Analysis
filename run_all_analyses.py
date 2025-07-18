#!/usr/bin/env python3
"""
Genomic Analysis Runner

This utility originally attempted to run a ``Pharmacy.py`` module which was not
included in the repository and therefore caused a ``FileNotFoundError``.  The
runner now executes only scripts that actually exist.  The available analyses
include ``Disease.py``, ``Comprehensive_Disease.py``, ``Ethnicity.py`` and
``Fitness_Athletics.py``.
"""

import subprocess
import sys

def run_analysis(script_name):
    """Run a Python analysis script and print its output"""
    print(f"\n{'='*60}")
    print(f"RUNNING {script_name.upper()} ANALYSIS")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run([sys.executable, script_name], 
                              capture_output=True, text=True, check=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        print(f"Error output: {e.stderr}")

if __name__ == "__main__":
    print("GENOMIC SEQUENCE ANALYSIS")
    print("Analyzing your genetic data...")
    
    # Run the available analyses
    run_analysis("Disease.py")
    run_analysis("Comprehensive_Disease.py")
    run_analysis("Ethnicity.py")
    run_analysis("Fitness_Athletics.py")
    
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print("Check the results above for your genetic insights!") 
