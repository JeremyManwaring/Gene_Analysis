"""End-to-end CLI test for one-command analyze output bundle."""

from __future__ import annotations

import json
from pathlib import Path

from Gene_Analysis.cli import main


def test_cli_analyze_writes_json_and_markdown(tmp_path: Path) -> None:
    genome_file = tmp_path / "sample_23andme.txt"
    genome_file.write_text(
        "\n".join(
            [
                "rs429358\t19\t45411941\tC\tC",
                "rs7903146\t10\t114758349\tC\tT",
                "rs4988235\t2\t136608646\tC\tC",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    snapshot_file = tmp_path / "snapshot.json"
    snapshot_file.write_text(
        json.dumps(
            {
                "metadata": {
                    "snapshot_id": "test-snapshot",
                    "created_at": "2026-04-01T00:00:00+00:00",
                    "stale_after_days": 30,
                },
                "markers": {},
                "trait_models": {},
                "gene_metadata": {},
            }
        ),
        encoding="utf-8",
    )

    json_out = tmp_path / "report.json"
    md_out = tmp_path / "report.md"
    rc = main(
        [
            "analyze",
            str(genome_file),
            "--no-auto-update",
            "--offline",
            "--include-legacy-catalog",
            "--cache-dir",
            str(tmp_path / "cache"),
            "--output-dir",
            str(tmp_path),
            "--snapshot-path",
            str(snapshot_file),
            "--json-output",
            str(json_out),
            "--markdown-output",
            str(md_out),
            "--top-traits",
            "5",
        ]
    )

    assert rc == 0
    assert json_out.exists()
    assert md_out.exists()

    report = json.loads(json_out.read_text(encoding="utf-8"))
    assert report["rows"], "Expected at least one SNP row in report"
    required_fields = {
        "trait",
        "rsid",
        "snp_type",
        "genotype",
        "risk_allele",
        "relative_probability",
        "error_rate",
        "confidence",
        "bias_note",
        "source",
    }
    assert required_fields.issubset(set(report["rows"][0].keys()))
