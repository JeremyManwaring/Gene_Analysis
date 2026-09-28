"""ClinVar fixtures validate allele/build matching, provenance and failure isolation."""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from Gene_Analysis.evidence import clinical


def make_vcf(path: Path, *, assembly: str = "GRCh37", release: str = "20260923") -> Path:
    records = [
        "1\t100\t111\tA\tG\t.\t.\tALLELEID=1;GENEINFO=ABC:1;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;CLNDN=Example_condition;CLNORIGIN=1",
        "1\t100\t112\tA\tT\t.\t.\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts",
        "1\t200\t113\tC\tT\t.\t.\tCLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Pathogenic(1)|Benign(1);CLNREVSTAT=criteria_provided,_conflicting_classifications",
        "1\t300\t114\tA\tAT\t.\t.\tCLNSIG=Pathogenic;CLNREVSTAT=practice_guideline",
        "X\t400\t115\tG\tT\t.\t.\tCLNSIG=Likely_pathogenic;CLNREVSTAT=criteria_provided,_single_submitter",
        "2\t500\t116\tC\tG\t.\t.\tCLNSIG=Uncertain_significance;CLNREVSTAT=no_assertion_criteria_provided",
    ]
    content = (f"##fileformat=VCFv4.2\n##fileDate={release}\n##reference={assembly}\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(records) + "\n")
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            handle.write(content)
    else:
        path.write_text(content)
    return path


def data(records, *, assembly="GRCh37", strand="plus"):
    return SimpleNamespace(assembly=assembly, strand=strand, get_dataframe=lambda: records)


def record(rsid, chromosome, position, genotype):
    return {"rsid": rsid, "chromosome": chromosome, "position": str(position), "genotype": genotype}


@pytest.fixture()
def manifest(tmp_path):
    source = make_vcf(tmp_path / "fixture.vcf.gz")
    return clinical.refresh_clinvar_database(source_path=source, cache_dir=tmp_path / "cache")


def test_import_and_matching_cover_exact_alleles_haploid_review_conflicts(manifest):
    genome = data([
        record("rsUnrelated", "1", 100, "AG"),
        record("rsAtWrongPosition", "1", 101, "AG"),
        record("rsConflict", "1", 200, "CT"),
        record("rsIndel", "1", 300, "AT"),
        record("rsHaploid", "23", 400, "T"),
        record("rsMismatch", "2", 500, "AG"),
    ])
    report = clinical.match_clinvar(genome, manifest_path=manifest["manifest_path"])
    assert report["status"] == "complete"
    assert {row["clinvar_variation_id"] for row in report["findings"]} == {"111", "113", "115"}
    first = report["findings"][0]
    assert first["review_stars"] == 3 and first["expert_reviewed"]
    assert first["needs_professional_review"]
    assert first["alternate_copies"] == 1
    assert first["absolute_risk"] is None and first["time_dependent_risk"] is None
    conflict = next(row for row in report["findings"] if row["clinvar_variation_id"] == "113")
    assert conflict["conflicting_classifications"] and not conflict["needs_professional_review"]
    assert report["summary"]["coordinates_without_clinvar"] == 1
    assert report["summary"]["allele_mismatches"] == 1
    assert report["summary"]["unsupported_variant_records"] == 1
    assert manifest["variant_count"] == 6 and manifest["snv_count"] == 5


def test_reference_is_not_variant_and_unknown_build_or_strand_withholds(manifest):
    genome = data([record("rs111", "1", 100, "AA")])
    report = clinical.match_clinvar(genome, manifest_path=manifest["manifest_path"])
    assert not report["findings"] and report["summary"]["reference_only_records"] == 1
    assert clinical.match_clinvar(data([], assembly="unknown"))["status"] == "withheld"
    assert clinical.match_clinvar(data([], strand="unknown"))["status"] == "withheld"
    assert clinical.match_clinvar(data([], strand="minus"))["status"] == "withheld"
    with pytest.raises(ValueError, match="assembly|build"):
        clinical.match_clinvar(genome, manifest_path=manifest["manifest_path"], assembly="GRCh38")
    with pytest.raises(ValueError, match="match requested"):
        clinical.match_clinvar(data([], assembly="GRCh38"), manifest_path=manifest["manifest_path"])


def test_manifest_roundtrip_hashes_and_local_source_is_not_live_verified(manifest):
    loaded = clinical.load_clinvar_manifest(manifest["manifest_path"])
    assert loaded == manifest
    assert loaded["fetched_at"] is None
    assert loaded["source_kind"] == "user_supplied_file"
    assert loaded["upstream_checksum_verified"] is False
    assert len(loaded["source_sha256"]) == len(loaded["database_sha256"]) == 64
    path = Path(manifest["manifest_path"])
    payload = json.loads(path.read_text())
    payload["source_release"] = "20260925"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="manifest checksum"):
        clinical.load_clinvar_manifest(path)


def test_corrupted_database_and_missing_pin_fail_without_network(manifest, tmp_path):
    with Path(manifest["database_path"]).open("ab") as handle:
        handle.write(b"changed")
    with pytest.raises(ValueError, match="database checksum"):
        clinical.load_clinvar_manifest(manifest["manifest_path"])
    with pytest.raises(FileNotFoundError, match="Pinned"):
        clinical.load_clinvar_manifest(tmp_path / "does-not-exist.json")


def test_failed_import_preserves_current_database_and_pinned_files(manifest, tmp_path):
    cache = Path(manifest["manifest_path"]).parent
    pointer = cache / "latest_clinvar_GRCh37.json"
    original = pointer.read_bytes()
    wrong = make_vcf(tmp_path / "wrong.vcf", assembly="GRCh38")
    with pytest.raises(ValueError, match="does not confirm"):
        clinical.refresh_clinvar_database(source_path=wrong, cache_dir=cache)
    assert pointer.read_bytes() == original
    assert clinical.load_clinvar_manifest(cache_dir=cache)["database_sha256"] == manifest["database_sha256"]
    empty = tmp_path / "empty.vcf"
    empty.write_text("##fileformat=VCFv4.2\n##fileDate=20260923\n##reference=GRCh37\n"
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with pytest.raises(ValueError, match="no variants"):
        clinical.refresh_clinvar_database(source_path=empty, cache_dir=cache)
    assert pointer.read_bytes() == original


def test_reimport_does_not_mutate_pinned_manifest(manifest, tmp_path):
    pinned_path = Path(manifest["manifest_path"])
    original = pinned_path.read_bytes()
    source = make_vcf(tmp_path / "new.vcf.gz", release="20260924")
    new = clinical.refresh_clinvar_database(source_path=source, cache_dir=pinned_path.parent)
    assert pinned_path.read_bytes() == original
    assert new["manifest_path"] != manifest["manifest_path"]
    assert clinical.load_clinvar_manifest(pinned_path)["source_release"] == "20260923"


def test_fresh_copy_of_old_release_remains_stale():
    assert clinical.is_clinvar_stale({"status": "complete", "source_release": "20010101"})
    assert clinical.is_clinvar_stale({"status": "partial", "source_release": "20260923"})
    with pytest.raises(ValueError):
        clinical.is_clinvar_stale(None, 0)


def test_checksum_failure_does_not_publish_download(monkeypatch, tmp_path):
    class Response(io.BytesIO):
        headers = {"Content-Length": "4"}
        def geturl(self):
            return "https://example.test/source"
    monkeypatch.setattr(clinical.urllib.request, "urlopen", lambda *args, **kwargs: Response(b"data"))
    path = tmp_path / "partial.gz"
    with pytest.raises(RuntimeError, match="checksum mismatch"):
        clinical._download("https://example.test/source", path, "0" * 32, attempts=1)
    assert not path.exists()


def test_review_stars_are_not_pathogenicity_probabilities():
    assert clinical.review_stars("practice_guideline") == 4
    assert clinical.review_stars("reviewed_by_expert_panel") == 3
    assert clinical.review_stars("criteria_provided,_multiple_submitters,_no_conflicts") == 2
    assert clinical.review_stars("criteria_provided,_conflicting_classifications") == 1
    assert clinical.review_stars("no_assertion_criteria_provided") == 0


def test_conflicting_calls_at_same_coordinate_are_not_interpreted(manifest):
    genome = data([record("rs1", "1", 100, "AG"), record("i123", "1", 100, "AA")])
    result = clinical.match_clinvar(genome, manifest_path=manifest["manifest_path"])
    assert result["summary"]["conflicting_coordinate_records"] == 2
    assert result["findings"] == []


def test_multiple_submitters_no_conflicts_is_not_a_conflict(tmp_path):
    source = make_vcf(tmp_path / "two-star.vcf")
    source.write_text(source.read_text().replace(
        "CLNSIG=Benign;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts",
        "CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts"))
    manifest = clinical.refresh_clinvar_database(source_path=source, cache_dir=tmp_path / "cache")
    report = clinical.match_clinvar(data([record("rs112", "1", 100, "AT")]),
                                   manifest_path=manifest["manifest_path"])
    finding = report["findings"][0]
    assert finding["review_stars"] == 2
    assert not finding["conflicting_classifications"]
    assert finding["needs_professional_review"]


def test_official_iso_file_date_is_accepted_without_changing_provenance(tmp_path):
    source = make_vcf(tmp_path / 'official-format.vcf', release='2026-09-23')
    manifest = clinical.refresh_clinvar_database(source_path=source, cache_dir=tmp_path / 'cache')
    assert manifest['source_release'] == '20260923'
    assert manifest['source_header']['fileDate'] == '2026-09-23'
    assert clinical.load_clinvar_manifest(manifest['manifest_path'], assembly=None)['assembly'] == 'GRCh37'


def test_complete_interrupted_download_is_reused_only_after_upstream_checksum(tmp_path, monkeypatch):
    import hashlib
    destination = tmp_path / 'source.gz'
    destination.write_bytes(b'complete public source')
    expected = hashlib.md5(destination.read_bytes()).hexdigest()
    requests = []
    class HeadResponse:
        headers = {'Last-Modified': 'Wed, 23 Sep 2026 00:00:00 GMT'}
        def __enter__(self): return self
        def __exit__(self, *args): return False
        def geturl(self): return 'https://example.invalid/clinvar.gz'
    def open_head(request, **kwargs):
        requests.append(request.get_method())
        assert request.get_method() == 'HEAD'
        return HeadResponse()
    monkeypatch.setattr(clinical.urllib.request, 'urlopen', open_head)
    provenance = clinical._download('https://example.invalid/clinvar.gz', destination, expected)
    assert requests == ['HEAD']
    assert provenance['upstream_checksum_verified'] is True
    assert provenance['source_md5'] == expected
