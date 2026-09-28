"""Input integrity and provenance, independent of personal genotypes."""
import gzip
from pathlib import Path

import pytest

from Gene_Analysis.utils import load_genome_data


def write_input(tmp_path: Path, content: str, name: str = "raw.txt") -> Path:
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    return path


def test_explicit_header_provenance_and_haploids(tmp_path):
    path = write_input(tmp_path, "# 23andMe\n# reference human genome build 37\n# plus strand\n"
                       "rs1\tchr1\t0001\tag\nrs2\tX\t2\tC\nrs3\tMT\t3\tT\nrs4\t1\t4\tA\n")
    data = load_genome_data(path)
    assert (data.assembly, data.strand, data.provider) == ("GRCh37", "plus", "23andMe")
    assert data.assembly_source == data.strand_source == "input_header"
    assert data.get_all_genotypes() == {"rs1": "AG", "rs2": "C", "rs3": "T"}
    assert data.stats.haploid_calls == 2
    assert data.stats.skipped_unsupported_genotypes == 1
    assert data.get_dataframe()[0]["position"] == "1"


def test_unknown_provenance_stays_unknown_and_explicit_override(tmp_path):
    path = write_input(tmp_path, "rs1 1 12 AG\n")
    data = load_genome_data(path)
    assert data.assembly == data.strand == "unknown"
    overridden = load_genome_data(path, assembly="hg38", strand="plus")
    assert overridden.assembly == "GRCh38"
    assert overridden.assembly_source == overridden.strand_source == "user_provided"


@pytest.mark.parametrize("header,options", [
    ("# build 37", {"assembly": "GRCh38"}),
    ("# plus strand", {"strand": "minus"}),
    ("# build 37\n# build 38", {}),
    ("# plus strand\n# minus strand", {}),
    ("# NCBI36", {}),
])
def test_conflicting_or_unsupported_provenance_is_rejected(tmp_path, header, options):
    path = write_input(tmp_path, header + "\nrs1\t1\t12\tAG\n")
    with pytest.raises(ValueError):
        load_genome_data(path, **options)


def test_ancestry_split_alleles_are_unphased_and_numeric_chromosomes(tmp_path):
    path = write_input(tmp_path, "\ufeff# AncestryDNA\n# reference genome version 37.1\n"
                       "rsid\tchromosome\tposition\tallele1\tallele2\n"
                       "rs1\t1\t12\tG\tA\nrs2\t23\t13\tC\tC\nrs3\t26\t14\tA\tA\n"
                       "rs4\t1\t15\t0\t0\nrs5\t1\t16\t-\tA\n")
    data = load_genome_data(path)
    assert data.assembly == "GRCh37"
    assert data.strand == "unknown"
    assert data.get_all_genotypes() == {"rs1": "AG", "rs2": "CC", "rs3": "AA"}
    assert [row["chromosome"] for row in data.get_dataframe()] == ["1", "X", "MT"]
    assert data.stats.skipped_header_rows == 1
    assert data.stats.skipped_no_calls == 2


def test_duplicate_conflicts_are_quarantined_not_last_write_wins(tmp_path):
    path = write_input(tmp_path, "rs1 1 1 AG\nrs1 1 1 GA\nrs2 1 2 CC\nrs2 1 2 CT\n"
                       "rs2 1 2 CC\nrs3 1 3 AA\nrs3 1 4 AA\nrs4 1 5 TT\n")
    data = load_genome_data(path)
    assert data.get_all_genotypes() == {"rs1": "AG", "rs4": "TT"}
    assert data.stats.parsed_rows == 2
    assert data.stats.duplicate_rows == 4
    assert data.stats.excluded_conflicting_markers == 2


def test_gzip_by_magic_and_utf8_bom(tmp_path):
    path = tmp_path / "renamed_raw.txt"
    path.write_bytes(gzip.compress("\ufeff# build 38\nrs1\t1\t10\tAA\n".encode("utf-8")))
    data = load_genome_data(path)
    assert data.get_genotype("rs1") == "AA"
    assert data.assembly == "GRCh38"


@pytest.mark.parametrize("content", ["", "# nothing here\n", "rs1 1 1 --\n", "rs1,1,1,AA\n",
                                      "bad 1 1 AA\n", "rs1 1 0 AA\n", "rs1 100 1 AA\n"])
def test_empty_malformed_and_no_usable_calls_fail_even_without_validate(tmp_path, content):
    path = write_input(tmp_path, content)
    with pytest.raises(ValueError, match="No valid"):
        load_genome_data(path, validate=False)


@pytest.mark.parametrize("content,name", [("##fileformat=VCFv4.3\n", "renamed.txt"),
    ("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", "renamed.txt"),
    ("rs1 1 1 AA\n", "sample.vcf"), ("rs1 1 1 AA\n", "sample.vcf.gz"),
    ("rs1 1 1 A G extra\n", "raw.txt")])
def test_unsupported_formats_are_not_silently_misparsed(tmp_path, content, name):
    path = write_input(tmp_path, content, name)
    with pytest.raises(ValueError, match="not supported|Unsupported"):
        load_genome_data(path)


def test_invalid_rows_and_indels_excluded_without_manufacturing_calls(tmp_path):
    path = write_input(tmp_path, "rs1 1 1 AA\nrs2 1 2 DD\nrs3 1 3 II\nrs4 1 4 NN\n"
                       "rs5 1 invalid AG\nrs6 1 6 A -\nrs7\t1\t7\t\tG\n")
    data = load_genome_data(path)
    assert data.get_all_genotypes() == {"rs1": "AA"}
    assert data.stats.skipped_no_calls == 3
    assert data.stats.skipped_unsupported_genotypes == 2
    assert data.stats.skipped_invalid_rows == 1


def test_returned_records_do_not_mutate_stored_calls(tmp_path):
    data = load_genome_data(write_input(tmp_path, "rs1 1 1 AA\n"))
    data.get_dataframe()[0]["genotype"] = "CC"
    assert data.get_dataframe()[0]["genotype"] == "AA"


@pytest.mark.parametrize("declaration", ["strand plus", "strand: plus", "strand = plus", "plus strand"])
def test_explicit_strand_declarations(tmp_path, declaration):
    data = load_genome_data(write_input(tmp_path, "# " + declaration + "\nrs1 1 1 AA\n"))
    assert data.strand == "plus"
    assert data.strand_source == "input_header"


def test_missing_second_allele_is_not_misread_as_a_haploid_call(tmp_path):
    data = load_genome_data(write_input(tmp_path, "rs1\t1\t1\tAA\nrs2\tX\t2\tA\t\n"))
    assert data.get_genotype("rs2") is None
    assert data.stats.skipped_no_calls == 1


def test_discordant_probes_at_the_same_locus_are_all_excluded(tmp_path):
    data = load_genome_data(write_input(tmp_path, "rs1 1 1 AA\ni20 1 1 AG\nrs3 1 3 TT\n"))
    assert data.get_all_genotypes() == {"rs3": "TT"}
    assert data.stats.excluded_conflicting_loci == 1
    assert data.stats.excluded_conflicting_markers == 2
