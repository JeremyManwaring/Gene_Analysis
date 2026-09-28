"""Output formatting and export utilities."""

import json
from datetime import datetime
from typing import Dict, Any, Optional, List, Tuple
from pathlib import Path


def format_console_output(results: Dict[str, Any], section_title: str = "RESULTS") -> str:
    """
    Format results for console output with consistent styling.
    
    Args:
        results: Dictionary containing analysis results
        section_title: Title for the section
    
    Returns:
        Formatted string for console output
    """
    output = []
    output.append(f"\n{'='*60}")
    output.append(section_title)
    output.append(f"{'='*60}\n")
    
    # Format different result types
    if isinstance(results, dict):
        for key, value in results.items():
            if isinstance(value, (dict, list)):
                output.append(f"{key}:")
                output.append(_format_nested(value, indent=2))
            else:
                output.append(f"{key}: {value}")
    
    return "\n".join(output)


def _format_nested(data: Any, indent: int = 0) -> str:
    """Format nested dictionaries and lists."""
    indent_str = " " * indent
    output = []
    
    if isinstance(data, dict):
        for key, value in data.items():
            if isinstance(value, (dict, list)):
                output.append(f"{indent_str}{key}:")
                output.append(_format_nested(value, indent + 2))
            else:
                output.append(f"{indent_str}{key}: {value}")
    elif isinstance(data, list):
        for item in data:
            if isinstance(item, (dict, list)):
                output.append(_format_nested(item, indent))
            else:
                output.append(f"{indent_str}- {item}")
    
    return "\n".join(output)


def export_to_json(
    all_results: Dict[str, Any],
    filename: Optional[str] = None,
    output_dir: str = "."
) -> str:
    """
    Export all analysis results to JSON file.
    
    Args:
        all_results: Dictionary containing all analysis results
        filename: Output filename (default: results_YYYYMMDD_HHMMSS.json)
        output_dir: Output directory (default: current directory)
    
    Returns:
        Path to created JSON file
    """
    # Generate filename if not provided
    if filename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"genetic_analysis_results_{timestamp}.json"
    
    # Ensure output directory exists
    output_path = Path(output_dir) / filename
    
    # Prepare export data
    export_data = {
        "metadata": {
            "analysis_date": datetime.now().isoformat(),
            "version": "1.0",
            "analysis_types": list(all_results.keys())
        },
        **all_results
    }
    
    # Write JSON file
    with open(output_path, 'w') as f:
        json.dump(export_data, f, indent=2, default=str)
    
    return str(output_path)


def create_summary_report(all_results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create a unified summary report from all analysis results.
    
    Args:
        all_results: Dictionary containing all analysis results
    
    Returns:
        Summary dictionary with key statistics
    """
    summary = {
        "total_analyses": len(all_results),
        "analysis_types": list(all_results.keys()),
        "summary_statistics": {}
    }
    
    # Extract summary statistics from each analysis
    for analysis_type, results in all_results.items():
        if isinstance(results, dict):
            # Look for common summary fields
            stats = {}
            
            if "summary" in results:
                stats = results["summary"]
            elif "statistics" in results:
                stats = results["statistics"]
            else:
                # Try to extract key metrics
                for key in ["snps_found", "total_snps", "coverage", "risk_level", "probability"]:
                    if key in results:
                        stats[key] = results[key]
            
            if stats:
                summary["summary_statistics"][analysis_type] = stats
    
    return summary


def format_snp_table(
    snp_details: List[Tuple[str, str, str, float, int, str]],
    headers: Optional[List[str]] = None
) -> str:
    """
    Format SNP details as a table.
    
    Args:
        snp_details: List of tuples (rsid, genotype, risk_allele, weight, count, description)
        headers: Optional custom headers
    
    Returns:
        Formatted table string
    """
    if not snp_details:
        return "No SNPs found."
    
    if headers is None:
        headers = ["SNP", "Genotype", "Risk", "Weight", "Count", "Description"]
    
    # Calculate column widths
    col_widths = [len(h) for h in headers]
    for row in snp_details:
        for i, cell in enumerate(row[:len(headers)]):
            col_widths[i] = max(col_widths[i], len(str(cell)))
    
    # Format header
    header_row = "  ".join(h.ljust(w) for h, w in zip(headers, col_widths))
    separator = "-" * len(header_row)
    
    # Format rows
    rows = [header_row, separator]
    for row in snp_details:
        row_str = "  ".join(str(cell).ljust(w) for cell, w in zip(row[:len(headers)], col_widths))
        rows.append(row_str)
    
    rows.append(separator)
    
    return "\n".join(rows)

