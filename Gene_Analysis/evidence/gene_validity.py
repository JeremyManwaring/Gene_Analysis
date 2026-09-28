"""Release-pinned ClinGen gene-disease validity context, not variant diagnoses."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
import re
import tempfile
import urllib.request
import uuid
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from .clinical import default_cache_dir, file_sha256

SOURCE_URL = "https://search.clinicalgenome.org/kb/gene-validity/download"
SCHEMA_VERSION = 1
SCOPE = ("All rows in the ClinGen Gene-Disease Validity export, including Limited, Disputed, "
         "Refuted and No Known Disease Relationship classifications. This is not all known genes "
         "or all clinically relevant diseases. Gene-disease validity does not establish the "
         "pathogenicity of an individual's variant or a probability of disease.")
COLUMNS = ["GENE SYMBOL", "GENE ID (HGNC)", "DISEASE LABEL", "DISEASE ID (MONDO)", "MOI", "SOP",
           "CLASSIFICATION", "ONLINE REPORT", "CLASSIFICATION DATE", "GCEP"]
KEYS = ["gene_symbol", "hgnc_id", "disease_label", "mondo_id", "mode_of_inheritance", "sop",
        "classification", "report_url", "classification_date", "expert_panel"]


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _digest(snapshot: dict[str, Any]) -> str:
    payload = {key: value for key, value in snapshot.items() if key != "snapshot_sha256"}
    payload["metadata"] = {key: value for key, value in payload["metadata"].items() if key != "snapshot_path"}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _atomic_json(path: Path, payload: dict, *, immutable: bool = False) -> None:
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if immutable:
            os.link(temporary, path)
        else:
            os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _parse_csv(raw: bytes) -> tuple[list[dict[str, str]], str]:
    reader = csv.reader(io.StringIO(raw.decode("utf-8-sig"), newline=""))
    found_header = False
    source_date = None
    records = []
    for line_number, row in enumerate(reader, 1):
        if not row or not any(value.strip() for value in row):
            continue
        first = row[0].strip()
        if first.startswith("FILE CREATED:"):
            value = first.partition(":")[2].strip()
            source_date = datetime.strptime(value, "%Y-%m-%d").replace(tzinfo=UTC).isoformat()
        if not found_header:
            if [value.strip() for value in row] == COLUMNS:
                found_header = True
            continue
        if all(not value.strip("+ ") for value in row):
            continue
        if len(row) != len(COLUMNS):
            raise ValueError(f"Invalid ClinGen CSV row width on line {line_number}.")
        record = dict(zip(KEYS, (value.strip() for value in row)))
        if (not record["gene_symbol"] or not re.fullmatch(r"HGNC:\d+", record["hgnc_id"])
                or not re.fullmatch(r"MONDO:\d+", record["mondo_id"])
                or not record["classification"]):
            raise ValueError(f"Missing or malformed ClinGen identifiers on line {line_number}.")
        records.append(record)
    if not found_header or not source_date or not records:
        raise ValueError("ClinGen source must contain FILE CREATED, the complete CSV header, and curations.")
    # Keep classifications verbatim, including refuting evidence and future categories.
    return records, source_date


def refresh_gene_validity(cache_dir: str | Path | None = None,
                          source_path: str | Path | None = None) -> dict[str, Any]:
    """Download/import the full export. Failed refreshes leave the previous pointer unchanged."""
    cache = Path(cache_dir).expanduser() if cache_dir is not None else default_cache_dir()
    cache.mkdir(parents=True, exist_ok=True)
    fetched_at = _now()
    source_metadata: dict[str, Any] = {"source_url": SOURCE_URL, "fetched_at": fetched_at,
                                       "source_kind": "local_import" if source_path is not None else "official_download"}
    if source_path is None:
        request = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "Gene-Analysis/1.0"})
        with urllib.request.urlopen(request, timeout=30) as response:
            raw = response.read(50_000_001)
            if len(raw) > 50_000_000:
                raise ValueError("ClinGen export exceeds the supported download size.")
            length = response.headers.get("Content-Length")
            if length is not None and len(raw) != int(length):
                raise ValueError("Incomplete ClinGen download (Content-Length mismatch).")
            source_metadata.update(source_url=response.geturl(), source_etag=response.headers.get("ETag"),
                                   source_last_modified=response.headers.get("Last-Modified"))
    else:
        raw = Path(source_path).read_bytes()
        source_metadata["imported_from"] = str(Path(source_path).resolve())
    records, source_date = _parse_csv(raw)
    source_hash = hashlib.sha256(raw).hexdigest()
    csv_path = cache / f"clingen-source-{source_hash}.csv"
    if csv_path.exists():
        if file_sha256(csv_path) != source_hash:
            raise ValueError("Stored ClinGen source checksum failed; refusing to overwrite evidence.")
    else:
        fd, temporary = tempfile.mkstemp(prefix=".clingen-source-", dir=cache)
        try:
            with os.fdopen(fd, "wb") as handle:
                handle.write(raw)
                handle.flush()
                os.fsync(handle.fileno())
            try:
                os.link(temporary, csv_path)
            except FileExistsError:
                if file_sha256(csv_path) != source_hash:
                    raise ValueError("Concurrent ClinGen source checksum mismatch.")
        finally:
            Path(temporary).unlink(missing_ok=True)
    classifications = Counter(record["classification"] for record in records)
    strong = [record for record in records if record["classification"] in {"Strong", "Definitive"}]
    metadata = {"schema_version": SCHEMA_VERSION, "snapshot_id": str(uuid.uuid4()),
                "created_at": fetched_at, **source_metadata, "source_created_at": source_date,
                "source_release": source_date[:10], "source_sha256": source_hash,
                "source_bytes": len(raw), "source_file": csv_path.name, "scope": SCOPE,
                "summary": {"genes": len({record["hgnc_id"] for record in records}),
                            "curations": len(records), "classification_counts": dict(classifications),
                            "strong_or_definitive_curations": len(strong),
                            "genes_with_strong_or_definitive_curation": len({r["hgnc_id"] for r in strong})}}
    snapshot = {"metadata": metadata, "curations": records}
    snapshot["snapshot_sha256"] = _digest(snapshot)
    path = cache / f"clingen-{source_date[:10]}-{metadata['snapshot_id']}.json"
    _atomic_json(path, snapshot, immutable=True)
    # Re-read and verify before publishing the newest usable snapshot pointer.
    result = load_gene_validity(path=path)
    _atomic_json(cache / "clingen-latest.json", {"snapshot_file": path.name, "file_sha256": file_sha256(path)})
    return result


def load_gene_validity(path: str | Path | None = None,
                       cache_dir: str | Path | None = None) -> dict[str, Any] | None:
    """Load a pinned snapshot and verify its hash; never download implicitly."""
    cache = Path(cache_dir).expanduser() if cache_dir is not None else default_cache_dir()
    candidate = Path(path).expanduser() if path is not None else cache / "clingen-latest.json"
    if not candidate.is_file():
        return None
    snapshot = json.loads(candidate.read_text(encoding="utf-8"))
    if "snapshot_file" in snapshot:
        filename = snapshot["snapshot_file"]
        if Path(filename).name != filename:
            raise ValueError("Invalid ClinGen snapshot pointer path.")
        target = candidate.parent / filename
        if not target.is_file() or file_sha256(target) != snapshot.get("file_sha256"):
            raise ValueError("ClinGen snapshot file checksum failed or snapshot is missing.")
        candidate = target
        snapshot = json.loads(candidate.read_text(encoding="utf-8"))
    metadata = snapshot.get("metadata", {})
    if metadata.get("schema_version") != SCHEMA_VERSION or not isinstance(snapshot.get("curations"), list):
        raise ValueError("Unsupported ClinGen snapshot schema.")
    if snapshot.get("snapshot_sha256") != _digest(snapshot):
        raise ValueError("ClinGen snapshot checksum failed; refusing modified evidence.")
    source_file = metadata.get("source_file")
    if source_file:
        if Path(source_file).name != source_file:
            raise ValueError("Invalid ClinGen source path.")
        source = candidate.parent / source_file
        if source.exists() and file_sha256(source) != metadata.get("source_sha256"):
            raise ValueError("ClinGen source checksum failed; refusing modified evidence.")
    metadata["snapshot_path"] = str(candidate.resolve())
    return snapshot


def is_gene_validity_stale(snapshot: dict[str, Any] | None, stale_after_days: int = 30) -> bool:
    """Age uses the source FILE CREATED date; copying an old source cannot refresh it."""
    if stale_after_days < 0:
        raise ValueError("stale_after_days must be nonnegative.")
    try:
        stamp = (snapshot or {}).get("metadata", {}).get("source_created_at")
        release = datetime.fromisoformat(stamp.replace("Z", "+00:00"))
        if release.tzinfo is None:
            release = release.replace(tzinfo=UTC)
        age = (datetime.now(UTC) - release).total_seconds() / 86400
        return age < -1 or age > stale_after_days
    except (AttributeError, TypeError, ValueError):
        return True


def attach_gene_validity(clinical: dict[str, Any], snapshot: dict[str, Any]) -> dict[str, Any]:
    """Annotate gene context. An exact MONDO match is labeled without asserting causality."""
    by_gene: dict[str, list[dict[str, str]]] = {}
    for record in snapshot.get("curations", []):
        by_gene.setdefault(record["gene_symbol"].upper(), []).append(record)
    for finding in clinical.get("findings", []):
        symbols = {part.partition(":")[0].strip().upper()
                   for part in re.split(r"[|,;]", finding.get("gene", "")) if part}
        mondo_ids = set(re.findall(r"MONDO:\d+", finding.get("condition_identifiers") or ""))
        contexts = []
        for symbol in sorted(symbols):
            for record in by_gene.get(symbol, []):
                exact = record["mondo_id"] in mondo_ids
                contexts.append({**record, "match_scope": "gene_and_exact_mondo_context" if exact else "gene_context_only",
                                 "disease_identifier_in_clinvar_aggregate": exact,
                                 "variant_pathogenicity_confirmed": False})
        finding["clingen_gene_validity"] = contexts
        finding["clingen_context_note"] = (
            "Gene-level evidence only. Even an exact disease identifier in an aggregate ClinVar record "
            "does not establish that this variant causes that disease or provide personal risk.")
    clinical["clingen_provenance"] = {**snapshot.get("metadata", {}), "snapshot_sha256": snapshot.get("snapshot_sha256")}
    return clinical
