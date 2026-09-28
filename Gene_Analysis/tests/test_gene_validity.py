"""ClinGen release provenance, integrity and non-diagnostic gene annotations."""
import csv
import json
from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

from Gene_Analysis.evidence.gene_validity import (
    COLUMNS, attach_gene_validity, is_gene_validity_stale, load_gene_validity, refresh_gene_validity,
)


def fixture_csv(path: Path, day: str = "2026-01-01") -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["CLINGEN GENE DISEASE VALIDITY CURATIONS"] + [""] * 9)
        writer.writerow(["FILE CREATED: " + day] + [""] * 9)
        writer.writerow(COLUMNS)
        writer.writerow(["++++"] * 10)
        for symbol, hgnc, mondo, classification in [
            ("TEST1", "HGNC:1", "MONDO:0000001", "Definitive"),
            ("TEST1", "HGNC:1", "MONDO:0000002", "Disputed"),
            ("TEST2", "HGNC:2", "MONDO:0000003", "Refuted"),
            ("TEST3", "HGNC:3", "MONDO:0000004", "Strong"),
        ]:
            writer.writerow([symbol, hgnc, "Example condition", mondo, "AD", "SOP10", classification,
                             "https://search.clinicalgenome.org/kb/gene-validity/example", "2025-01-01", "Example panel"])
    return path


def test_import_preserves_all_classes_and_reloads_by_pinned_or_latest(tmp_path):
    source = fixture_csv(tmp_path / "input.csv")
    cache = tmp_path / "cache"
    snapshot = refresh_gene_validity(cache_dir=cache, source_path=source)
    summary = snapshot["metadata"]["summary"]
    assert summary["curations"] == 4 and summary["genes"] == 3
    assert summary["strong_or_definitive_curations"] == 2
    assert summary["classification_counts"]["Disputed"] == 1
    assert summary["classification_counts"]["Refuted"] == 1
    assert snapshot["metadata"]["source_kind"] == "local_import"
    pinned = Path(snapshot["metadata"]["snapshot_path"])
    assert load_gene_validity(path=pinned) == load_gene_validity(cache_dir=cache) == snapshot
    assert load_gene_validity(path=tmp_path / "missing.json") is None


def test_source_age_not_import_age_controls_staleness(tmp_path):
    old_day = (datetime.now(UTC) - timedelta(days=60)).date().isoformat()
    old = refresh_gene_validity(cache_dir=tmp_path / "cache", source_path=fixture_csv(tmp_path / "old.csv", old_day))
    assert old["metadata"]["fetched_at"].startswith(datetime.now(UTC).date().isoformat())
    assert is_gene_validity_stale(old, 30)
    assert not is_gene_validity_stale(old, 90)
    assert is_gene_validity_stale(None)
    assert is_gene_validity_stale({"metadata": {"created_at": datetime.now(UTC).isoformat()}})


def test_failed_refresh_keeps_previous_snapshot_pointer(tmp_path):
    cache = tmp_path / "cache"
    good = refresh_gene_validity(cache_dir=cache, source_path=fixture_csv(tmp_path / "good.csv"))
    pointer_before = (cache / "clingen-latest.json").read_bytes()
    bad = tmp_path / "bad.csv"
    bad.write_text("<html>Service temporarily unavailable</html>")
    with pytest.raises(ValueError):
        refresh_gene_validity(cache_dir=cache, source_path=bad)
    assert (cache / "clingen-latest.json").read_bytes() == pointer_before
    assert load_gene_validity(cache_dir=cache) == good


def test_modified_snapshot_rejected_when_loading_directly(tmp_path):
    snapshot = refresh_gene_validity(cache_dir=tmp_path / "cache", source_path=fixture_csv(tmp_path / "source.csv"))
    path = Path(snapshot["metadata"]["snapshot_path"])
    payload = json.loads(path.read_text())
    payload["curations"][0]["classification"] = "Refuted"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="checksum"):
        load_gene_validity(path=path)


def test_modified_source_rejected(tmp_path):
    cache = tmp_path / "cache"
    snapshot = refresh_gene_validity(cache_dir=cache, source_path=fixture_csv(tmp_path / "source.csv"))
    (cache / snapshot["metadata"]["source_file"]).write_text("changed")
    with pytest.raises(ValueError, match="source checksum"):
        load_gene_validity(cache_dir=cache)


def test_gene_match_alone_never_confirms_variant_or_disease(tmp_path):
    snapshot = refresh_gene_validity(cache_dir=tmp_path / "cache", source_path=fixture_csv(tmp_path / "source.csv"))
    clinical = {"findings": [{"gene": "TEST1:1|TEST2:2", "condition_identifiers": "MONDO:0000001,MedGen:123"}]}
    result = attach_gene_validity(clinical, snapshot)
    annotations = result["findings"][0]["clingen_gene_validity"]
    assert len(annotations) == 3
    assert annotations[0]["match_scope"] == "gene_and_exact_mondo_context"
    assert annotations[1]["match_scope"] == "gene_context_only"
    assert annotations[1]["classification"] == "Disputed"
    assert annotations[2]["classification"] == "Refuted"
    assert not any(row["variant_pathogenicity_confirmed"] for row in annotations)
