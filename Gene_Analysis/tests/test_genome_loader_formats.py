"""Tests for 23andMe/raw genome parser format support."""

from __future__ import annotations

from pathlib import Path

from Gene_Analysis.utils import load_genome_data


def test_loads_4_column_genotype_format(tmp_path: Path) -> None:
    genome_file = tmp_path / "genome4.txt"
    genome_file.write_text(
        "\n".join(
            [
                "# comment",
                "rs1\t1\t1000\tAG",
                "rs2\t1\t1001\tCC",
                "rs3\t1\t1002\t--",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    data = load_genome_data(str(genome_file))
    assert data.get_genotype("rs1") == "AG"
    assert data.get_genotype("rs2") == "CC"
    assert data.get_genotype("rs3") is None
    assert data.stats.inferred_format == "4col_genotype"
    assert data.stats.parsed_rows == 2


def test_loads_5_column_phased_format(tmp_path: Path) -> None:
    genome_file = tmp_path / "genome5.txt"
    genome_file.write_text(
        "\n".join(
            [
                "rs10\t1\t2000\tA\tG",
                "rs11\t1\t2001\tC\tC",
                "rs12\t1\t2002\t-\tA",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    data = load_genome_data(str(genome_file))
    assert data.get_genotype("rs10") == "AG"
    assert data.get_genotype("rs11") == "CC"
    assert data.get_genotype("rs12") is None
    assert data.stats.inferred_format == "5col_allele1_allele2"
    assert data.stats.parsed_rows == 2


def test_large_file_smoke_parse(tmp_path: Path) -> None:
    genome_file = tmp_path / "large_genome.txt"
    lines = [f"rs{i}\t1\t{i}\tA\tG" for i in range(1, 25001)]
    genome_file.write_text("\n".join(lines) + "\n", encoding="utf-8")

    data = load_genome_data(str(genome_file))
    assert data.stats.parsed_rows == 25000
    assert data.get_genotype("rs24999") == "AG"
