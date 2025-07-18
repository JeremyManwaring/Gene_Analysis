#!/usr/bin/env python3
"""
Genomic Analysis Runner

This script is a simple wrapper that sequentially executes the various
analysis modules.  Earlier versions called ``Pharmacy.py`` which no longer
exists.  The runner now expects consolidated modules to be available and
will call them when present.
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
    
    # Invoke each analysis module. These files should exist in the same
    # directory as this runner script.  They may be implemented later, so
    # missing files will simply raise an error when executed.
    run_analysis("disease_analysis.py")
    run_analysis("ethnicity_analysis.py")
    run_analysis("longevity_analysis.py")
    
    print(f"\n{'='*60}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print("Check the results above for your genetic insights!") 