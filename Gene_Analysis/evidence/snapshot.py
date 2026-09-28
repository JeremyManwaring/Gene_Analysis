"""Evidence snapshot update/load utilities for external genetics sources."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import math
import re
import uuid
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from ..analysis.catalog import MarkerDefinition, build_marker_catalog, marker_rsids
from .http import HttpClient
from .clinical import default_cache_dir


GWAS_SNP_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms/{rsid}"
GWAS_ASSOC_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"
DBSNP_SUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
CLINVAR_TABLES_URL = "https://clinicaltables.nlm.nih.gov/api/variants/v4/search"
GNOMAD_GRAPHQL = "https://gnomad.broadinstitute.org/api"
PGS_TRAIT_ALL_URL = "https://www.pgscatalog.org/rest/trait/all"
PGS_SCORE_SEARCH_URL = "https://www.pgscatalog.org/rest/score/search"
PGS_PERFORMANCE_SEARCH_URL = "https://www.pgscatalog.org/rest/performance/search"
OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
CACHE_DIR = default_cache_dir() / "supplementary"
LATEST_PATH = CACHE_DIR / "latest_snapshot.json"


def _now_utc() -> datetime:
    return datetime.now(UTC)


def _safe_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if not text:
        return None
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    if not match:
        return None
    try:
        return float(match.group(0))
    except ValueError:
        return None


def _parse_isoish(value: Optional[str]) -> Optional[datetime]:
    if not value:
        return None
    text = value.strip()
    if not text:
        return None

    for candidate in [text.replace("Z", "+00:00"), text]:
        try:
            dt = datetime.fromisoformat(candidate)
            if dt.tzinfo is None:
                return dt.replace(tzinfo=UTC)
            return dt.astimezone(UTC)
        except ValueError:
            pass

    for fmt in ["%Y/%m/%d %H:%M", "%Y-%m-%d", "%Y/%m/%d"]:
        try:
            dt = datetime.strptime(text, fmt)
            return dt.replace(tzinfo=UTC)
        except ValueError:
            continue
    return None


def _clean_tokens(text: str) -> List[str]:
    cleaned = re.sub(r"[^a-z0-9]+", " ", text.lower())
    stopwords = {
        "risk",
        "disease",
        "late",
        "onset",
        "and",
        "or",
        "with",
        "for",
        "the",
        "of",
        "childhood",
        "region",
        "proxy",
        "deficiency",
        "metabolism",
    }
    return [tok for tok in cleaned.split() if tok and tok not in stopwords]


def _trait_overlap_score(trait_name: str, candidate: Dict[str, Any]) -> float:
    source_tokens = set(_clean_tokens(trait_name))
    if not source_tokens:
        return 0.0

    target_bits = [candidate.get("label", "")]
    target_bits.extend(candidate.get("trait_synonyms", []))
    target_tokens: set[str] = set()
    for bit in target_bits:
        target_tokens.update(_clean_tokens(bit))

    if not target_tokens:
        return 0.0

    inter = len(source_tokens.intersection(target_tokens))
    if inter == 0:
        return 0.0
    return inter / math.sqrt(len(source_tokens) * len(target_tokens))


def _payload_hash(payload: Dict[str, Any]) -> str:
    copy = {**payload, "metadata": dict(payload.get("metadata", {}))}
    copy["metadata"].pop("payload_sha256", None)
    return hashlib.sha256(json.dumps(copy, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _json_dump(path: Path, payload: Dict[str, Any], *, immutable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".snapshot-", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.flush()
            os.fsync(handle.fileno())
        if immutable:
            os.link(temporary, path)
        else:
            os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def load_snapshot(path: Optional[str] = None) -> Optional[Dict[str, Any]]:
    snapshot_path = Path(path) if path else LATEST_PATH
    if not snapshot_path.exists():
        if path:
            raise FileNotFoundError(f"Pinned evidence snapshot not found: {snapshot_path}")
        return None
    with snapshot_path.open("r", encoding="utf-8") as handle:
        snapshot = json.load(handle)
    expected = snapshot.get("metadata", {}).get("payload_sha256")
    if expected and expected != _payload_hash(snapshot):
        raise ValueError("Evidence snapshot checksum failed; refusing modified evidence.")
    return snapshot


def snapshot_age_days(snapshot: Optional[Dict[str, Any]]) -> Optional[float]:
    if not snapshot:
        return None
    created_at = _parse_isoish(snapshot.get("metadata", {}).get("created_at"))
    if not created_at:
        return None
    age = _now_utc() - created_at
    return age.total_seconds() / 86400.0


def is_snapshot_stale(snapshot: Optional[Dict[str, Any]], stale_after_days: int) -> bool:
    if stale_after_days < 1:
        raise ValueError("stale_after_days must be at least one.")
    if not snapshot or snapshot.get("metadata", {}).get("refresh_status") != "complete":
        return True
    age_days = snapshot_age_days(snapshot)
    if age_days is None or age_days < 0:
        return True
    return age_days > stale_after_days


def _fetch_gwas_snp(client: HttpClient, rsid: str) -> Dict[str, Any]:
    response = client.get_json(GWAS_SNP_URL.format(rsid=rsid))
    if response.get("_error"):
        return {"error": response["_error"]}
    return {
        "rsid": response.get("rs_id", rsid),
        "functional_class": response.get("functional_class"),
        "most_severe_consequence": response.get("most_severe_consequence"),
        "mapped_genes": response.get("mapped_genes", []),
        "alleles": response.get("alleles"),
        "last_update_date": response.get("last_update_date"),
        "locations": response.get("locations", []),
        "source": "GWAS Catalog",
        "source_url": GWAS_SNP_URL.format(rsid=rsid),
    }


def _fetch_gwas_associations(client: HttpClient, rsid: str, max_records: int = 30) -> List[Dict[str, Any]]:
    response = client.get_json(
        GWAS_ASSOC_URL,
        params={"rs_id": rsid, "page": 0, "size": max_records},
    )
    if response.get("_error"):
        return []

    rows = response.get("_embedded", {}).get("associations", [])
    parsed: List[Dict[str, Any]] = []
    for row in rows:
        parsed.append(
            {
                "association_id": row.get("association_id"),
                "p_value": _safe_float(row.get("p_value")),
                "or_value": _safe_float(row.get("or_value") or row.get("or_per_copy_num")),
                "beta": _safe_float(row.get("beta")),
                "ci_lower": _safe_float(row.get("ci_lower")),
                "ci_upper": _safe_float(row.get("ci_upper")),
                "risk_frequency": _safe_float(row.get("risk_frequency")),
                "pubmed_id": str(row.get("pubmed_id")) if row.get("pubmed_id") else None,
                "first_author": row.get("first_author"),
                "reported_trait": row.get("reported_trait", []),
                "efo_traits": row.get("efo_traits", []),
                "mapped_genes": row.get("mapped_genes", []),
                "snp_allele": row.get("snp_allele", []),
                "source": "GWAS Catalog",
            }
        )

    parsed.sort(key=lambda item: item.get("p_value") if item.get("p_value") is not None else 10.0)
    return parsed


def _fetch_dbsnp_summary(client: HttpClient, rsid: str) -> Dict[str, Any]:
    rs_num = rsid.lower().replace("rs", "", 1)
    if not rs_num.isdigit():
        return {"error": "non-rsid marker"}

    response = client.get_json(
        DBSNP_SUMMARY_URL,
        params={"db": "snp", "id": rs_num, "retmode": "json"},
    )
    if response.get("_error"):
        return {"error": response["_error"]}

    summary = response.get("result", {}).get(rs_num, {})
    if not summary:
        return {"error": "no dbSNP summary"}

    mafs = summary.get("global_mafs", []) or []
    gnomad_mafs = [m for m in mafs if "gnomad" in str(m.get("study", "")).lower()]
    return {
        "snp_id": summary.get("snp_id"),
        "genes": [gene.get("name") for gene in summary.get("genes", []) if gene.get("name")],
        "fxn_class": summary.get("fxn_class"),
        "spdi": summary.get("spdi"),
        "clinical_significance": summary.get("clinical_significance"),
        "global_mafs": mafs[:20],
        "gnomad_mafs": gnomad_mafs,
        "createdate": summary.get("createdate"),
        "updatedate": summary.get("updatedate"),
        "source": "dbSNP",
        "source_url": f"{DBSNP_SUMMARY_URL}?db=snp&id={rs_num}&retmode=json",
    }


def _fetch_clinvar_tables(client: HttpClient, rsid: str) -> Dict[str, Any]:
    response = client.get_json(CLINVAR_TABLES_URL, params={"terms": rsid})
    if isinstance(response, dict) and response.get("_error"):
        return {"error": response["_error"]}

    if not isinstance(response, list) or len(response) < 4:
        return {"error": "invalid Clinical Tables response shape"}

    count = int(response[0]) if str(response[0]).isdigit() else 0
    top_hits = response[3][:5] if len(response) > 3 else []
    parsed_hits = []
    for hit in top_hits:
        if isinstance(hit, list) and len(hit) >= 2:
            parsed_hits.append({"id": hit[0], "name": hit[1]})
    return {
        "count": count,
        "top_hits": parsed_hits,
        "source": "ClinVar Clinical Tables",
        "clinical_interpretation_allowed": False,
        "match_scope": "rsID search only; neither genotype nor alternate allele verified",
        "source_url": f"{CLINVAR_TABLES_URL}?terms={rsid}",
    }


def _fetch_gnomad(client: HttpClient, rsid: str) -> Dict[str, Any]:
    search_payload = {
        "query": (
            "query ($q: String!) { "
            "variant_search(query: $q, dataset: gnomad_r4) { variant_id }"
            "}"
        ),
        "variables": {"q": rsid},
    }
    search_resp = client.post_json(GNOMAD_GRAPHQL, search_payload)
    if search_resp.get("_error"):
        return {"error": search_resp["_error"]}

    hits = search_resp.get("data", {}).get("variant_search", [])
    if not hits:
        return {"error": "variant_not_found"}

    variant_id = hits[0].get("variant_id")
    if not variant_id:
        return {"error": "variant_not_found"}

    detail_payload = {
        "query": (
            "query ($variantId: String!) { "
            "variant(variantId: $variantId, dataset: gnomad_r4) { "
            "variantId "
            "exome { ac an af populations { id ac an } } "
            "genome { ac an af populations { id ac an } } "
            "}"
            "}"
        ),
        "variables": {"variantId": variant_id},
    }
    detail_resp = client.post_json(GNOMAD_GRAPHQL, detail_payload)
    if detail_resp.get("_error"):
        return {"error": detail_resp["_error"], "variant_id": variant_id}

    variant = detail_resp.get("data", {}).get("variant")
    if not variant:
        return {"error": "variant_details_missing", "variant_id": variant_id}

    exome_block = variant.get("exome") or {}
    genome_block = variant.get("genome") or {}

    def _top_pops(block: Optional[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not block:
            return []
        rows = block.get("populations", [])
        scored = []
        for row in rows:
            ac = _safe_float(row.get("ac")) or 0.0
            an = _safe_float(row.get("an")) or 0.0
            af = ac / an if an > 0 else 0.0
            scored.append({"id": row.get("id"), "ac": ac, "an": an, "af": af})
        scored.sort(key=lambda item: item["af"], reverse=True)
        return scored[:10]

    return {
        "variant_id": variant.get("variantId"),
        "exome_af": _safe_float(exome_block.get("af")),
        "genome_af": _safe_float(genome_block.get("af")),
        "exome_top_populations": _top_pops(exome_block if isinstance(exome_block, dict) else None),
        "genome_top_populations": _top_pops(genome_block if isinstance(genome_block, dict) else None),
        "source": "gnomAD",
        "source_url": GNOMAD_GRAPHQL,
    }


def _fetch_pgs_trait_index(client: HttpClient) -> List[Dict[str, Any]]:
    traits: List[Dict[str, Any]] = []
    page = 0
    size = 250
    while True:
        resp = client.get_json(PGS_TRAIT_ALL_URL, params={"limit": size, "offset": page * size})
        if resp.get("_error"):
            break
        rows = resp.get("results", [])
        if not rows:
            break
        traits.extend(rows)
        if not resp.get("next"):
            break
        page += 1
    return traits


def _pick_best_pgs_trait(trait_name: str, trait_index: List[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    if not trait_index:
        return None

    best = None
    best_score = 0.0
    for candidate in trait_index:
        score = _trait_overlap_score(trait_name, candidate)
        if score > best_score:
            best_score = score
            best = candidate

    if best and best_score >= 0.25:
        return best
    return None


def _parse_date_key(value: Optional[str]) -> datetime:
    dt = _parse_isoish(value)
    if dt is None:
        return datetime(1970, 1, 1, tzinfo=UTC)
    return dt


def _fetch_pgs_model_for_trait(client: HttpClient, trait_name: str, trait_index: List[Dict[str, Any]]) -> Dict[str, Any]:
    trait_match = _pick_best_pgs_trait(trait_name, trait_index)
    if not trait_match:
        return {"error": "trait_not_matched"}

    trait_id = trait_match.get("id")
    score_resp = client.get_json(PGS_SCORE_SEARCH_URL, params={"trait_id": trait_id, "limit": 8})
    if score_resp.get("_error"):
        return {"error": score_resp["_error"], "trait_id": trait_id}

    scores = score_resp.get("results", [])
    if not scores:
        return {"error": "no_scores", "trait_id": trait_id}

    best_score = max(scores, key=lambda item: _parse_date_key(item.get("date_release")))
    pgs_id = best_score.get("id")

    perf_resp = client.get_json(PGS_PERFORMANCE_SEARCH_URL, params={"pgs_id": pgs_id, "limit": 10})
    if perf_resp.get("_error"):
        performances: List[Dict[str, Any]] = []
    else:
        performances = perf_resp.get("results", [])

    best_perf = None
    if performances:
        best_perf = max(
            performances,
            key=lambda item: _parse_date_key(item.get("publication", {}).get("date_publication")),
        )

    return {
        "trait_name": trait_name,
        "trait_match": {
            "id": trait_match.get("id"),
            "label": trait_match.get("label"),
            "synonyms": trait_match.get("trait_synonyms", []),
        },
        "score": {
            "id": best_score.get("id"),
            "name": best_score.get("name"),
            "date_release": best_score.get("date_release"),
            "variants_number": best_score.get("variants_number"),
            "method_name": best_score.get("method_name"),
            "trait_reported": best_score.get("trait_reported"),
            "ancestry_distribution": best_score.get("ancestry_distribution", {}),
            "publication": best_score.get("publication", {}),
            "source_url": f"{PGS_SCORE_SEARCH_URL}?trait_id={trait_id}",
        },
        "performance": {
            "id": best_perf.get("id") if best_perf else None,
            "metrics": best_perf.get("performance_metrics", {}) if best_perf else {},
            "samples": (best_perf.get("sampleset", {}) or {}).get("samples", [])[:5] if best_perf else [],
            "publication": best_perf.get("publication", {}) if best_perf else {},
            "source_url": f"{PGS_PERFORMANCE_SEARCH_URL}?pgs_id={pgs_id}",
        },
        "source": "PGS Catalog",
    }


def _fetch_open_targets_gene(client: HttpClient, gene_symbol: str) -> Dict[str, Any]:
    search_payload = {
        "query": (
            "query ($query: String!) { "
            "search(queryString: $query, entityNames:[\"target\"], page:{index:0,size:1}) { "
            "total hits { id entity object { ... on Target { approvedSymbol approvedName } } }"
            "}"
            "}"
        ),
        "variables": {"query": gene_symbol},
    }
    search_resp = client.post_json(OT_GRAPHQL, search_payload)
    if search_resp.get("_error"):
        return {"error": search_resp["_error"]}

    hits = (search_resp.get("data", {}).get("search", {}) or {}).get("hits", [])
    if not hits:
        return {"error": "target_not_found"}

    hit = hits[0]
    ensembl_id = hit.get("id")
    target_payload = {
        "query": (
            "query ($ensemblId: String!) { "
            "target(ensemblId: $ensemblId) { "
            "id approvedSymbol approvedName "
            "associatedDiseases(page:{index:0,size:1}) { count rows { disease { id name } score } }"
            "}"
            "}"
        ),
        "variables": {"ensemblId": ensembl_id},
    }
    target_resp = client.post_json(OT_GRAPHQL, target_payload)
    if target_resp.get("_error"):
        return {
            "ensembl_id": ensembl_id,
            "approved_symbol": (hit.get("object") or {}).get("approvedSymbol"),
            "approved_name": (hit.get("object") or {}).get("approvedName"),
            "error": target_resp["_error"],
        }

    target = target_resp.get("data", {}).get("target", {})
    assoc = target.get("associatedDiseases", {})
    top_row = (assoc.get("rows") or [None])[0]
    return {
        "ensembl_id": target.get("id"),
        "approved_symbol": target.get("approvedSymbol"),
        "approved_name": target.get("approvedName"),
        "associated_disease_count": assoc.get("count"),
        "top_associated_disease": top_row,
        "source": "Open Targets",
        "source_url": OT_GRAPHQL,
    }


def _best_marker_timestamp(marker_payload: Dict[str, Any]) -> Optional[str]:
    candidates = [
        marker_payload.get("gwas_snp", {}).get("last_update_date"),
        marker_payload.get("dbsnp", {}).get("updatedate"),
    ]
    parsed = [dt for dt in (_parse_isoish(x) for x in candidates) if dt is not None]
    if not parsed:
        return None
    return max(parsed).isoformat()


def update_evidence_snapshot(
    markers: Optional[List[MarkerDefinition]] = None,
    stale_after_days: int = 30,
    output_path: Optional[str] = None,
    show_progress: bool = False,
) -> Dict[str, Any]:
    """Fetch and cache latest evidence snapshot for marker rsids/traits."""
    markers = build_marker_catalog() if markers is None else markers
    rsids = marker_rsids(markers)
    client = HttpClient(timeout_seconds=8, retries=1, backoff_seconds=0.4)
    gnomad_client = HttpClient(timeout_seconds=6, retries=0, backoff_seconds=0.2)
    ot_client = HttpClient(timeout_seconds=8, retries=0, backoff_seconds=0.2)
    rsid_categories: Dict[str, set[str]] = {}
    for marker in markers:
        rsid_categories.setdefault(marker.rsid, set()).add(marker.category)

    snapshot: Dict[str, Any] = {
        "metadata": {
            "snapshot_id": str(uuid.uuid4()),
            "schema_version": 2,
            "scope": "Supplementary public source metadata; not clinical allele interpretation or risk calibration",
            "created_at": _now_utc().isoformat(),
            "stale_after_days": stale_after_days,
            "rsid_count": len(rsids),
            "traits_count": len({m.trait for m in markers if m.category == "disease"}),
            "sources": [
                "GWAS Catalog",
                "PGS Catalog",
                "dbSNP",
                "ClinVar Clinical Tables",
                "gnomAD",
                "Open Targets",
            ],
        },
        "markers": {},
        "trait_models": {},
        "gene_metadata": {},
    }

    if show_progress:
        print(f"Refreshing marker evidence for {len(rsids)} unique rsIDs...")

    for index, rsid in enumerate(rsids, start=1):
        if show_progress and (index == 1 or index % 10 == 0 or index == len(rsids)):
            print(f"  - Marker {index}/{len(rsids)}: {rsid}")
        try:
            gwas_snp = _fetch_gwas_snp(client, rsid)
            gwas_associations = _fetch_gwas_associations(client, rsid, max_records=20)
            dbsnp = _fetch_dbsnp_summary(client, rsid)
            clinvar = _fetch_clinvar_tables(client, rsid)
            marker_categories = rsid_categories.get(rsid, set())
            if marker_categories.intersection({"disease", "carrier_screening", "pharmacogenomics"}):
                gnomad = _fetch_gnomad(gnomad_client, rsid)
            else:
                gnomad = {"skipped": True, "reason": "category_not_selected"}

            snapshot["markers"][rsid] = {
                "gwas_snp": gwas_snp,
                "gwas_associations": gwas_associations,
                "dbsnp": dbsnp,
                "clinvar": clinvar,
                "gnomad": gnomad,
            }
            snapshot["markers"][rsid]["best_last_updated"] = _best_marker_timestamp(snapshot["markers"][rsid])
        except Exception as exc:
            snapshot["markers"][rsid] = {
                "error": f"{type(exc).__name__}: {exc}",
                "gwas_snp": {},
                "gwas_associations": [],
                "dbsnp": {},
                "clinvar": {},
                "gnomad": {},
                "best_last_updated": None,
            }

    disease_traits = sorted({m.trait for m in markers if m.category == "disease"})
    pgs_trait_index = _fetch_pgs_trait_index(client)
    if show_progress:
        print(f"Refreshing PRS model metadata for {len(disease_traits)} disease traits...")
    for index, trait in enumerate(disease_traits, start=1):
        if show_progress and (index == 1 or index % 8 == 0 or index == len(disease_traits)):
            print(f"  - Trait {index}/{len(disease_traits)}: {trait}")
        try:
            snapshot["trait_models"][trait] = _fetch_pgs_model_for_trait(client, trait, pgs_trait_index)
        except Exception as exc:
            snapshot["trait_models"][trait] = {"error": f"{type(exc).__name__}: {exc}"}

    all_genes = sorted({m.gene for m in markers if m.gene and m.gene.isalpha() and len(m.gene) <= 10})
    if show_progress:
        print(f"Refreshing Open Targets metadata for {len(all_genes)} genes...")
    for index, gene in enumerate(all_genes, start=1):
        if show_progress and (index == 1 or index % 12 == 0 or index == len(all_genes)):
            print(f"  - Gene {index}/{len(all_genes)}: {gene}")
        try:
            snapshot["gene_metadata"][gene] = _fetch_open_targets_gene(ot_client, gene)
        except Exception as exc:
            snapshot["gene_metadata"][gene] = {"error": f"{type(exc).__name__}: {exc}"}

    snapshot["metadata"]["updated_at"] = _now_utc().isoformat()
    requests = client.requests + gnomad_client.requests + ot_client.requests
    failures = [request for request in requests if request["status"] != "success"]
    successes = len(requests) - len(failures)
    # A transport success is insufficient when parsers detected unusable data.
    benign_absence = {"trait_not_matched", "no_scores", "target_not_found", "variant_not_found",
                     "variant_details_missing", "no dbSNP summary", "non-rsid marker"}
    parser_errors = []
    def find_errors(value: Any, location: str = "") -> None:
        if isinstance(value, dict):
            if value.get("error") and value["error"] not in benign_absence:
                parser_errors.append({"location": location, "error": value["error"]})
            for key, nested in value.items():
                find_errors(nested, f"{location}/{key}")
        elif isinstance(value, list):
            for index, nested in enumerate(value):
                find_errors(nested, f"{location}/{index}")
    find_errors(snapshot)
    complete = bool(requests) and not failures and not parser_errors
    snapshot["metadata"].update(
        refresh_status="complete" if complete else ("partial" if successes else "failed"),
        requests_succeeded=successes, requests_failed=len(failures),
        parser_errors=parser_errors, request_provenance=requests,
        freshness_note="Retrieval time is not a clinical validation or a source-release date.",
    )
    target_path = Path(output_path) if output_path else CACHE_DIR / (
        f"snapshot_{snapshot['metadata']['snapshot_id']}.json")
    snapshot["metadata"]["path"] = str(target_path.resolve())
    snapshot["metadata"]["payload_sha256"] = _payload_hash(snapshot)
    _json_dump(target_path, snapshot, immutable=True)
    # Failed/partial refreshes remain auditable but cannot displace a complete cache.
    if complete:
        _json_dump(LATEST_PATH, snapshot)
    return snapshot
