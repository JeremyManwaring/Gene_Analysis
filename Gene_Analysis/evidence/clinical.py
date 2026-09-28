"""Local, release-pinned ClinVar annotation; never calculates disease probabilities.

ClinVar's VCF is a *subset* of ClinVar (mapped simple alleles under 10 kb).
Only exact, plus-strand single-nucleotide alleles are matched to raw genotypes.
No genotype data is sent over the network.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import sqlite3
import tempfile
import time
import urllib.request
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.parse import unquote

SCHEMA_VERSION = 1
CLINVAR_BASE_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_{assembly}/clinvar.vcf.gz"
CLINVAR_DOCS = "https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/"
BASES = frozenset("ACGT")
SCOPE = (
    "ClinVar VCF: precisely mapped simple alleles under 10 kb. Raw-genotype matching "
    "supports single-nucleotide variants only. Structural variants, repeats, indels, "
    "haplotypes, phase, penetrance and time-dependent disease probabilities are not evaluated."
)


def default_cache_dir() -> Path:
    """Writable user cache; installed package directories are never modified."""
    configured = os.environ.get("GENE_ANALYSIS_CACHE_DIR")
    return Path(configured).expanduser() if configured else Path.home() / ".cache" / "gene-analysis"


def _now() -> str:
    return datetime.now(UTC).isoformat()


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_hash(value: dict[str, Any]) -> str:
    payload = dict(value)
    payload.pop("manifest_sha256", None)
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _atomic_json(path: Path, value: dict[str, Any], *, immutable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        if immutable:
            os.link(temporary, path)  # fails rather than overwriting a pinned artifact
        else:
            os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _check_assembly(assembly: str) -> str:
    mapping = {"37": "GRCh37", "hg19": "GRCh37", "GRCh37": "GRCh37",
               "38": "GRCh38", "hg38": "GRCh38", "GRCh38": "GRCh38"}
    if assembly not in mapping:
        raise ValueError("ClinVar requires a known genome assembly: GRCh37 or GRCh38.")
    return mapping[assembly]


def _download(url: str, destination: Path, expected_md5: str, *, attempts: int = 3) -> dict[str, Any]:
    """Bounded retries, full download integrity check; never publishes a partial file."""
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            digest = hashlib.md5()  # NCBI publishes MD5; SHA256 is also recorded locally
            downloaded = 0
            # Resume a partial public download; the upstream checksum validates all bytes.
            if destination.exists():
                with destination.open("rb") as previous:
                    for block in iter(lambda: previous.read(1024 * 1024), b""):
                        digest.update(block)
                        downloaded += len(block)
            if downloaded and digest.hexdigest() == expected_md5:
                # Reuse a complete public download after an interrupted build.
                request = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "Gene-Analysis/1.0"})
                with urllib.request.urlopen(request, timeout=30) as response:
                    return {"source_url": response.geturl(), "fetched_at": _now(),
                            "source_last_modified": response.headers.get("Last-Modified"),
                            "source_etag": response.headers.get("ETag"), "source_bytes": downloaded,
                            "source_md5": digest.hexdigest(), "source_sha256": file_sha256(destination),
                            "upstream_checksum_verified": True}
            started = time.monotonic()
            headers = {"User-Agent": "Gene-Analysis/1.0"}
            if downloaded:
                headers["Range"] = f"bytes={downloaded}-"
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=30) as response:
                if downloaded and response.status != 206:
                    downloaded = 0
                    digest = hashlib.md5()
                resumed_bytes = downloaded
                length = response.headers.get("Content-Length")
                provenance = {
                    "source_url": response.geturl(), "fetched_at": _now(),
                    "source_last_modified": response.headers.get("Last-Modified"),
                    "source_etag": response.headers.get("ETag"),
                }
                with destination.open("ab" if resumed_bytes else "wb") as handle:
                    while block := response.read(1024 * 1024):
                        downloaded += len(block)
                        if downloaded > 2_000_000_000 or time.monotonic() - started > 1800:
                            raise ValueError("ClinVar download exceeded the size/time limit.")
                        digest.update(block)
                        handle.write(block)
                    handle.flush()
                    os.fsync(handle.fileno())
            if length and downloaded - resumed_bytes != int(length):
                raise ValueError("Incomplete ClinVar download (Content-Length mismatch).")
            if digest.hexdigest() != expected_md5:
                raise ValueError("ClinVar checksum mismatch; the source may have changed during download.")
            return {**provenance, "source_bytes": downloaded, "source_md5": digest.hexdigest(),
                    "source_sha256": file_sha256(destination), "upstream_checksum_verified": True}
        except (OSError, ValueError) as exc:
            last_error = exc
            if isinstance(exc, ValueError) or getattr(exc, "code", None) == 416:
                destination.unlink(missing_ok=True)
            if attempt + 1 < attempts:
                time.sleep(attempt + 1)
    raise RuntimeError(f"ClinVar download failed after {attempts} attempts: {last_error}") from last_error


def _read_checksum(url: str) -> str:
    request = urllib.request.Request(url, headers={"User-Agent": "Gene-Analysis/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        text = response.read(8192).decode("ascii")
    match = re.match(r"^([a-fA-F0-9]{32})\s", text)
    if not match:
        raise ValueError("NCBI did not provide a valid ClinVar MD5 checksum.")
    return match.group(1).lower()


def _chromosome(value: str) -> str:
    value = value.removeprefix("chr").upper()
    return {"23": "X", "24": "Y", "M": "MT", "26": "MT"}.get(value, value)


def _info_fields(text: str) -> dict[str, str]:
    return dict(token.split("=", 1) for token in text.split(";") if "=" in token)


def _vcf_text(value: str) -> str:
    # VCF field escapes are percent encoded, and ClinVar uses underscores for spaces.
    return unquote(value).replace("_", " ")


def _create_database(source: Path, temporary_db: Path, assembly: str) -> dict[str, Any]:
    source_header: dict[str, str] = {}
    count = 0
    snvs = 0
    opener = gzip.open if source.suffix == ".gz" else open
    with sqlite3.connect(temporary_db) as conn:
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
        conn.execute("CREATE TABLE variants (chromosome TEXT, position INTEGER, reference TEXT, "
                     "alternate TEXT, variation_id TEXT, details TEXT)")
        batch = []
        header_seen = False
        with opener(source, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("##"):
                    key, _, value = line[2:].strip().partition("=")
                    if key in {"fileDate", "reference", "source", "fileformat"}:
                        source_header[key] = value
                    continue
                if line.startswith("#CHROM\t"):
                    header_seen = True
                    reference = source_header.get("reference", "")
                    if assembly not in reference:
                        raise ValueError(f"VCF reference {reference!r} does not confirm {assembly}.")
                    if not re.fullmatch(r"(?:\d{8}|\d{4}-\d{2}-\d{2})", source_header.get("fileDate", "")):
                        raise ValueError("ClinVar VCF needs a fileDate release identifier.")
                    continue
                if line.startswith("#") or not line.strip():
                    continue
                if not header_seen:
                    raise ValueError("Missing VCF column header.")
                columns = line.rstrip("\r\n").split("\t")
                if len(columns) < 8:
                    raise ValueError(f"Malformed ClinVar VCF row {count + 1}.")
                chrom, position, variant_id, ref, alt = columns[:5]
                if not position.isdigit() or int(position) < 1 or not variant_id.isdigit():
                    raise ValueError("Invalid ClinVar coordinate or Variation ID.")
                info = _info_fields(columns[7])
                details = {key: info[key] for key in (
                    "ALLELEID", "GENEINFO", "CLNSIG", "CLNREVSTAT", "CLNDN", "CLNDISDB",
                    "CLNSIGCONF", "CLNORIGIN", "CLNVC", "RS", "CLNSIGINCL", "CLNDNINCL",
                    "ONC", "ONCREVSTAT", "SCI", "SCIREVSTAT",
                ) if key in info}
                batch.append((_chromosome(chrom), int(position), ref, alt, variant_id,
                              json.dumps(details, separators=(",", ":"))))
                count += 1
                snvs += ref in BASES and alt in BASES
                if len(batch) >= 10000:
                    conn.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?)", batch)
                    batch.clear()
            if batch:
                conn.executemany("INSERT INTO variants VALUES (?,?,?,?,?,?)", batch)
        if count == 0:
            raise ValueError("ClinVar VCF contains no variants; refusing empty database.")
        source_release = source_header["fileDate"].replace("-", "")
        # Parse real calendar dates, not merely an eight digit string.
        datetime.strptime(source_release, "%Y%m%d")
        conn.execute("CREATE INDEX variant_coordinates ON variants(chromosome, position)")
        metadata = {"assembly": assembly, "source_release": source_release,
                    "schema_version": SCHEMA_VERSION, "variant_count": count, "snv_count": snvs}
        conn.executemany("INSERT INTO metadata VALUES (?,?)",
                         [(key, json.dumps(value)) for key, value in metadata.items()])
        conn.commit()
        if conn.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise ValueError("ClinVar SQLite integrity verification failed.")
    return {**metadata, "source_header": source_header}


def refresh_clinvar_database(
    assembly: str = "GRCh37", cache_dir: str | Path | None = None,
    source_path: str | Path | None = None, show_progress: bool = False,
    *, resume_path: str | Path | None = None,
) -> dict[str, Any]:
    """Download/import a complete ClinVar VCF and atomically publish a pinned database.

    Local source imports are offline and marked as user supplied, not live-verified.
    The current pointer is updated only after checksum, parse and SQLite validation.
    """
    assembly = _check_assembly(assembly)
    cache = Path(cache_dir) if cache_dir else default_cache_dir()
    cache.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".clinvar-build-", dir=cache) as temporary:
        temporary_dir = Path(temporary)
        if source_path is None:
            url = CLINVAR_BASE_URL.format(assembly=assembly)
            if show_progress:
                print(f"Downloading public ClinVar {assembly} VCF (about 200 MB); genotypes stay local...")
            expected_md5 = _read_checksum(url + ".md5")
            downloads = cache / "downloads"
            downloads.mkdir(exist_ok=True)
            source = Path(resume_path) if resume_path else downloads / f"clinvar_{assembly}.vcf.gz"
            provenance = _download(url, source, expected_md5)
            provenance["source_kind"] = "official_ncbi_download"
        else:
            source = Path(source_path).expanduser().resolve()
            if not source.is_file():
                raise FileNotFoundError(f"ClinVar VCF not found: {source}")
            provenance = {"source_url": None, "source_kind": "user_supplied_file",
                          "source_sha256": file_sha256(source), "source_bytes": source.stat().st_size,
                          "fetched_at": None, "imported_at": _now(), "upstream_checksum_verified": False}
        if show_progress:
            print("Building and verifying the local ClinVar database...")
        temporary_db = temporary_dir / "clinvar.sqlite3"
        metadata = _create_database(source, temporary_db, assembly)
        if file_sha256(source) != provenance["source_sha256"]:
            raise ValueError("ClinVar source changed during import; database was not published.")
        db_hash = file_sha256(temporary_db)
        stem = f"clinvar_{assembly}_{metadata['source_release']}_{db_hash[:16]}"
        database_path = cache / f"{stem}.sqlite3"
        if database_path.exists():
            if file_sha256(database_path) != db_hash:
                raise ValueError("Existing pinned ClinVar database has changed; refusing to replace it.")
        else:
            os.link(temporary_db, database_path)
        manifest = {
            **metadata, **provenance, "created_at": _now(), "status": "complete",
            "database_path": database_path.name, "database_sha256": db_hash,
            "scope": SCOPE, "source_documentation": CLINVAR_DOCS,
        }
        manifest["manifest_sha256"] = _canonical_hash(manifest)
        manifest_path = cache / f"{stem}_{manifest['manifest_sha256'][:12]}.json"
        _atomic_json(manifest_path, manifest, immutable=True)
        _atomic_json(cache / f"latest_clinvar_{assembly}.json", {
            "manifest_path": manifest_path.name, "manifest_sha256": manifest["manifest_sha256"]})
    # Return portable data plus resolved locations; the saved manifest stays immutable.
    return {**manifest, "database_path": str(database_path.resolve()),
            "manifest_path": str(manifest_path.resolve())}


def load_clinvar_manifest(
    path: str | Path | None = None, assembly: str | None = "GRCh37",
    *, verify_database: bool = True, cache_dir: str | Path | None = None,
) -> dict[str, Any] | None:
    """Resolve an exact manifest (or latest pointer), verify hashes, and open locally."""
    if assembly is not None:
        assembly = _check_assembly(assembly)
    elif path is None:
        raise ValueError("A genome assembly is required when selecting a cached ClinVar database.")
    if path is None:
        cache = Path(cache_dir) if cache_dir else default_cache_dir()
        candidate = cache / f"latest_clinvar_{assembly}.json"
        if not candidate.exists():
            return None
    else:
        candidate = Path(path).expanduser()
        if not candidate.is_file():
            raise FileNotFoundError(f"Pinned ClinVar manifest not found: {candidate}")
    payload = json.loads(candidate.read_text(encoding="utf-8"))
    pointer_hash = None
    if "manifest_path" in payload and "database_path" not in payload:
        pointer_hash = payload.get("manifest_sha256")
        candidate = candidate.parent / payload["manifest_path"]
        payload = json.loads(candidate.read_text(encoding="utf-8"))
    if payload.get("status") != "complete" or payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("ClinVar manifest is incomplete or uses an unsupported schema.")
    expected = payload.get("manifest_sha256")
    if expected != _canonical_hash(payload) or (pointer_hash and expected != pointer_hash):
        raise ValueError("ClinVar manifest checksum failed; refusing modified evidence.")
    source_assembly = _check_assembly(payload.get("assembly"))
    if assembly is not None and source_assembly != assembly:
        raise ValueError(f"ClinVar {payload.get('assembly')} does not match requested {assembly}.")
    database_path = candidate.parent / payload["database_path"]
    if not database_path.is_file():
        raise FileNotFoundError(f"ClinVar database not found: {database_path}")
    if verify_database and file_sha256(database_path) != payload["database_sha256"]:
        raise ValueError("ClinVar database checksum failed; refusing modified evidence.")
    return {**payload, "manifest_path": str(candidate.resolve()),
            "database_path": str(database_path.resolve())}


def is_clinvar_stale(manifest: dict[str, Any] | None, stale_after_days: int = 30) -> bool:
    """Age is based on the source release, never a recent copy/import timestamp."""
    if stale_after_days < 1:
        raise ValueError("stale_after_days must be at least one.")
    if not manifest or manifest.get("status") != "complete":
        return True
    try:
        release = datetime.strptime(manifest["source_release"], "%Y%m%d").replace(tzinfo=UTC)
    except (KeyError, TypeError, ValueError):
        return True
    age = (datetime.now(UTC) - release).total_seconds() / 86400
    return age < -1 or age > stale_after_days


def review_stars(review_status: str) -> int:
    normalized = review_status.replace("_", " ").lower()
    if normalized == "practice guideline":
        return 4
    if normalized == "reviewed by expert panel":
        return 3
    if normalized == "criteria provided, multiple submitters, no conflicts":
        return 2
    if normalized in {"criteria provided, single submitter", "criteria provided, conflicting classifications",
                       "criteria provided, conflicting interpretations"}:
        return 1
    return 0


def _finding(row: sqlite3.Row, record: dict[str, str]) -> dict[str, Any]:
    info = json.loads(row["details"])
    significance = _vcf_text(info.get("CLNSIG", "not provided"))
    review = _vcf_text(info.get("CLNREVSTAT", "not provided"))
    stars = review_stars(review)
    conflict = (bool(info.get("CLNSIGCONF")) or "conflict" in significance.lower()
                or "conflicting" in review.lower())
    classes = {part.strip().lower() for part in re.split(r"[/|,]", significance)}
    pathogenic = bool(classes & {"pathogenic", "likely pathogenic"})
    return {
        "rsid": record.get("rsid"), "chromosome": row["chromosome"], "position": row["position"],
        "reference": row["reference"], "alternate": row["alternate"], "genotype": record["genotype"],
        "alternate_copies": record["genotype"].count(row["alternate"]),
        "clinvar_variation_id": row["variation_id"], "clinvar_allele_id": info.get("ALLELEID"),
        "gene": info.get("GENEINFO", ""), "germline_classification": significance,
        "review_status": review, "review_stars": stars,
        "expert_reviewed": stars >= 3, "conflicting_classifications": conflict,
        "needs_professional_review": pathogenic and not conflict and stars >= 2,
        "conditions": _vcf_text(info.get("CLNDN", "not provided")),
        "condition_identifiers": info.get("CLNDISDB"), "allele_origin_codes": info.get("CLNORIGIN"),
        "classification_scope": "variant aggregate; individual condition assertions require source review",
        "match_method": "exact assembly, chromosome, position and observed plus-strand alternate SNV",
        "source_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{row['variation_id']}/",
        "absolute_risk": None, "time_dependent_risk": None,
    }


def match_clinvar(
    genome_data: Any, manifest_path: str | Path | None = None, assembly: str | None = None,
    *, manifest: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Annotate observed alternate SNVs, requiring known build and plus strand.

    This does not interpret inheritance, diagnose disease, assert penetrance, or
    interpret absent/unmeasured alleles. A correct allele match is not assay validation.
    """
    input_assembly = getattr(genome_data, "assembly", "unknown")
    input_strand = getattr(genome_data, "strand", "unknown")
    result: dict[str, Any] = {
        "status": "withheld", "scope": SCOPE, "assembly": input_assembly, "findings": [],
        "limitations": [
            "ClinVar submissions are not independently verified by NIH; review status varies.",
            "Consumer raw genotypes require clinical confirmation before medical decisions.",
            "No match does not rule out any condition; most disease-causing variation may be unmeasured.",
            "A variant classification is not a personal disease probability or an inheritance assessment.",
        ],
    }
    if input_assembly not in {"GRCh37", "GRCh38"}:
        result["reason"] = "Unknown input genome build; supply the build confirmed by the data provider."
        return result
    if assembly is not None and _check_assembly(assembly) != input_assembly:
        raise ValueError("Requested ClinVar assembly differs from the input genome build.")
    if input_strand != "plus":
        result["reason"] = "Input plus-strand orientation is not confirmed; clinical matching withheld."
        return result
    selected_manifest = manifest_path or (manifest.get("manifest_path") if manifest else None)
    manifest = load_clinvar_manifest(selected_manifest, assembly=input_assembly)
    if manifest is None:
        result["reason"] = "No local ClinVar database. Run update-evidence for the input genome build."
        return result
    if manifest.get("assembly") != input_assembly:
        raise ValueError("ClinVar database assembly differs from the input genome build.")
    result["database"] = {key: manifest.get(key) for key in (
        "source_release", "source_sha256", "database_sha256", "manifest_sha256", "fetched_at",
        "source_url", "source_kind", "manifest_path", "variant_count", "snv_count",
    )}
    result["database"]["stale"] = is_clinvar_stale(manifest)
    count = Counter({key: 0 for key in (
        "input_records", "unsupported_genotypes", "invalid_coordinates",
        "unsupported_variant_records", "allele_mismatches", "coordinates_without_clinvar",
        "records_with_match", "reference_only_records", "conflicting_coordinate_records",
    )})
    findings = []
    seen = set()
    database = Path(manifest["database_path"])
    records = genome_data.get_dataframe()
    coordinate_genotypes: dict[tuple[str, str], str] = {}
    conflicting_coordinates = set()
    for record in records:
        key = (_chromosome(str(record.get("chromosome", ""))), str(record.get("position", "")))
        genotype = "".join(sorted(record.get("genotype", "")))
        if key in coordinate_genotypes and coordinate_genotypes[key] != genotype:
            conflicting_coordinates.add(key)
        coordinate_genotypes[key] = genotype
    with sqlite3.connect(database.as_uri() + "?mode=ro", uri=True) as conn:
        conn.row_factory = sqlite3.Row
        for record in records:
            count["input_records"] += 1
            genotype = record.get("genotype", "")
            if len(genotype) not in {1, 2} or any(base not in BASES for base in genotype):
                count["unsupported_genotypes"] += 1
                continue
            chromosome = _chromosome(str(record.get("chromosome", "")))
            position = str(record.get("position", ""))
            if (chromosome, position) in conflicting_coordinates:
                count["conflicting_coordinate_records"] += 1
                continue
            if not position.isdigit() or int(position) <= 0:
                count["invalid_coordinates"] += 1
                continue
            rows = conn.execute("SELECT * FROM variants WHERE chromosome=? AND position=?",
                                (chromosome, int(position))).fetchall()
            if not rows:
                count["coordinates_without_clinvar"] += 1
            matched_record = False
            reference_only = False
            for row in rows:
                if row["reference"] not in BASES or row["alternate"] not in BASES:
                    count["unsupported_variant_records"] += 1
                    continue
                if row["alternate"] not in genotype:
                    reference_only |= set(genotype) == {row["reference"]}
                    continue
                if not set(genotype).issubset({row["reference"], row["alternate"]}):
                    count["allele_mismatches"] += 1
                    continue
                identity = (chromosome, row["position"], row["reference"], row["alternate"], row["variation_id"])
                if identity in seen:
                    continue
                seen.add(identity)
                findings.append(_finding(row, record))
                matched_record = True
            count["records_with_match"] += matched_record
            count["reference_only_records"] += reference_only and not matched_record
    findings.sort(key=lambda item: (not item["needs_professional_review"], -item["review_stars"],
                                    item["chromosome"], item["position"], item["clinvar_variation_id"]))
    count["matched_variants"] = len(findings)
    count["expert_reviewed_matches"] = sum(item["expert_reviewed"] for item in findings)
    count["conflicting_matches"] = sum(item["conflicting_classifications"] for item in findings)
    count["professional_review_candidates"] = sum(item["needs_professional_review"] for item in findings)
    result.update(status="complete", summary=dict(count), findings=findings)
    return result
