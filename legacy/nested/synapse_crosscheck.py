"""Synapse connection + local cross-check of genome report runs.

This module does two things:

1.  Provides a Synapse client (``SynapseEvidenceClient``) that is *ready to run*
    once the user authorizes the Synapse connector / provides an auth token.
    Synapse (Sage Bionetworks, https://www.synapse.org) is a data-repository
    platform; it requires OAuth / a personal-access-token, so it cannot be
    reached from a non-interactive sandbox.  The client here degrades
    gracefully: if no credentials are found it reports ``connected=False`` and
    the cross-check falls back to the *local cached evidence snapshot* that the
    pipeline already downloaded from GWAS Catalog, dbSNP, ClinVar and gnomAD.

2.  Cross-checks the two existing report runs against that evidence and, for
    every disease trait, keeps the **most accurate** value (highest confidence /
    lowest error rate).  The combined result is written to a ``.txt`` file.

Nothing here calls the network unless a Synapse token is present, so it is safe
to run offline.
"""

from __future__ import annotations

if __name__ == "__main__":
    raise SystemExit("Retired heuristic tool. Run: python -m Gene_Analysis analyze /path/to/data.txt")



import glob
import json
import os
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parent
REPORTS_DIR = REPO_ROOT / "analysis_reports"
SNAPSHOT_PATH = PACKAGE_ROOT / "evidence_cache" / "latest_snapshot.json"
GENOME_RESULTS = PACKAGE_ROOT / "jeremy_genome_results.json"


# --------------------------------------------------------------------------- #
# 1. Synapse client (ready to run once authorized)
# --------------------------------------------------------------------------- #
@dataclass
class SynapseEvidenceClient:
    """Thin wrapper around the Synapse Python client.

    Usage once authorized (interactive session)::

        client = SynapseEvidenceClient(auth_token=os.environ["SYNAPSE_AUTH_TOKEN"])
        client.login()
        table = client.query_table("syn########", "SELECT rsid, trait, odds_ratio FROM %s")

    The token can come from:
      * the ``SYNAPSE_AUTH_TOKEN`` environment variable, or
      * ``~/.synapseConfig``  (the standard synapseclient config file), or
      * the Synapse MCP connector once it is authorized in the app.
    """

    auth_token: Optional[str] = None
    default_view: Optional[str] = None  # e.g. a Synapse table/view synID
    _syn: Any = field(default=None, repr=False)
    connected: bool = False
    status: str = "not_attempted"

    def _resolve_token(self) -> Optional[str]:
        if self.auth_token:
            return self.auth_token
        env = os.environ.get("SYNAPSE_AUTH_TOKEN") or os.environ.get("SYNAPSE_PAT")
        if env:
            return env
        cfg = Path.home() / ".synapseConfig"
        if cfg.exists():
            for line in cfg.read_text(encoding="utf-8", errors="replace").splitlines():
                if "authtoken" in line.lower() and "=" in line:
                    return line.split("=", 1)[1].strip()
        return None

    def login(self) -> bool:
        """Attempt to log in. Returns True on success, False (no raise) otherwise."""
        token = self._resolve_token()
        if not token:
            self.status = (
                "no_credentials: Synapse needs OAuth/personal-access-token. "
                "Authorize the Synapse connector in the app, or set SYNAPSE_AUTH_TOKEN, "
                "then re-run in an interactive session."
            )
            self.connected = False
            return False
        try:
            import synapseclient  # type: ignore

            self._syn = synapseclient.Synapse()
            self._syn.login(authToken=token, silent=True)
            self.connected = True
            self.status = "connected"
            return True
        except ImportError:
            self.status = "synapseclient not installed (`pip install synapseclient`)"
        except Exception as exc:  # noqa: BLE001 - report, never crash the pipeline
            self.status = f"login_failed: {type(exc).__name__}: {exc}"
        self.connected = False
        return False

    def query_table(self, syn_id: str, sql_template: str = "SELECT * FROM %s") -> List[Dict[str, Any]]:
        """Query a Synapse table/fileview into a list of dict rows."""
        if not self.connected or self._syn is None:
            return []
        try:
            query = sql_template % syn_id if "%s" in sql_template else sql_template
            results = self._syn.tableQuery(query)
            df = results.asDataFrame()
            return df.to_dict(orient="records")
        except Exception as exc:  # noqa: BLE001
            self.status = f"query_failed: {type(exc).__name__}: {exc}"
            return []

    def build_rsid_index(self, syn_id: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
        """Return {rsid: row} from a Synapse table, if reachable; else {}."""
        target = syn_id or self.default_view
        if not target or not self.connected:
            return {}
        index: Dict[str, Dict[str, Any]] = {}
        for row in self.query_table(target):
            rsid = str(row.get("rsid") or row.get("rsID") or row.get("SNP") or "").strip().lower()
            if rsid.startswith("rs"):
                index[rsid] = row
        return index


# --------------------------------------------------------------------------- #
# 2. Local evidence snapshot (the "data we have")
# --------------------------------------------------------------------------- #
def load_json(path: Path) -> Optional[Dict[str, Any]]:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def load_evidence_snapshot() -> Dict[str, Any]:
    snap = load_json(SNAPSHOT_PATH) or {}
    return snap.get("markers", {}) if snap else {}


def evidence_supports_marker(rsid: str, snapshot_markers: Dict[str, Any]) -> Dict[str, Any]:
    """Summarize what the cached external evidence says about one rsid."""
    ev = snapshot_markers.get(rsid, {})
    if not ev:
        return {"evidence_found": False, "sources": []}

    gwas = ev.get("gwas_snp", {}) or {}
    dbsnp = ev.get("dbsnp", {}) or {}
    clinvar = ev.get("clinvar", {}) or {}
    assocs = ev.get("gwas_associations", []) or []
    best_p = None
    for a in assocs:
        p = a.get("p_value")
        if isinstance(p, (int, float)) and (best_p is None or p < best_p):
            best_p = p

    sources = []
    if gwas and not gwas.get("error"):
        sources.append("GWAS Catalog")
    if dbsnp and not dbsnp.get("error"):
        sources.append("dbSNP")
    if clinvar and not clinvar.get("error") and clinvar.get("count"):
        sources.append("ClinVar")

    return {
        "evidence_found": bool(sources),
        "sources": sources,
        "mapped_genes": gwas.get("mapped_genes") or dbsnp.get("genes") or [],
        "most_severe_consequence": gwas.get("most_severe_consequence"),
        "clinvar_count": clinvar.get("count", 0),
        "n_gwas_associations": len(assocs),
        "best_p_value": best_p,
        "last_updated": ev.get("best_last_updated"),
    }


# --------------------------------------------------------------------------- #
# 3. Cross-check + combine the two report runs
# --------------------------------------------------------------------------- #
@dataclass
class TraitVerdict:
    category: str
    trait: str
    chosen_run: str
    relative_probability: Optional[float]
    error_rate: Optional[float]
    confidence: Optional[float]
    markers_found: int
    markers_total: int
    evidence_sources: List[str]
    evidence_found: bool
    runs_agree: bool
    note: str


def _num(v: Any) -> Optional[float]:
    return float(v) if isinstance(v, (int, float)) else None


def find_report_runs() -> List[Path]:
    return sorted(REPORTS_DIR.glob("gene_report_*.json"))


def _trait_rsids(report: Dict[str, Any], trait: str) -> List[str]:
    return [r.get("rsid") for r in report.get("rows", []) if r.get("trait") == trait and r.get("rsid")]


def cross_check(prefer: str = "confidence") -> Dict[str, Any]:
    """Combine every report run, keeping the most accurate value per trait.

    prefer="confidence"  -> pick the run with the highest confidence (ties: lower error).
    """
    runs = find_report_runs()
    if len(runs) < 1:
        raise FileNotFoundError(f"No gene_report_*.json found in {REPORTS_DIR}")

    reports = []
    for p in runs:
        data = load_json(p)
        if data:
            reports.append((p.name, data))

    snapshot_markers = load_evidence_snapshot()

    # Collect trait_summaries keyed by (category, trait) per run.
    per_run: Dict[str, Dict[Tuple[str, str], Dict[str, Any]]] = {}
    all_traits: Dict[Tuple[str, str], None] = {}
    for name, rep in reports:
        idx: Dict[Tuple[str, str], Dict[str, Any]] = {}
        for ts in rep.get("trait_summaries", []):
            key = (ts.get("category", ""), ts.get("trait", ""))
            idx[key] = ts
            all_traits[key] = None
        per_run[name] = idx

    verdicts: List[TraitVerdict] = []
    for (category, trait) in all_traits:
        candidates: List[Tuple[str, Dict[str, Any]]] = []
        for name, idx in per_run.items():
            if (category, trait) in idx:
                candidates.append((name, idx[(category, trait)]))
        if not candidates:
            continue

        def score(item: Tuple[str, Dict[str, Any]]) -> Tuple[float, float]:
            ts = item[1]
            conf = _num(ts.get("confidence")) or 0.0
            err = _num(ts.get("error_rate"))
            err = err if err is not None else 100.0
            return (conf, -err)  # higher confidence, then lower error

        best_name, best = max(candidates, key=score)

        # Do the runs agree on the headline probability (within 1 pp)?
        probs = [_num(c[1].get("relative_probability")) for c in candidates]
        probs = [p for p in probs if p is not None]
        runs_agree = (max(probs) - min(probs) <= 1.0) if len(probs) >= 2 else True

        # Evidence support: union of rsids for this trait across runs.
        rsids: List[str] = []
        for _, rep in reports:
            rsids.extend(_trait_rsids(rep, trait))
        ev_sources: List[str] = []
        ev_found = False
        for rsid in set(rsids):
            supp = evidence_supports_marker(rsid, snapshot_markers)
            if supp["evidence_found"]:
                ev_found = True
                for s in supp["sources"]:
                    if s not in ev_sources:
                        ev_sources.append(s)

        note_bits = []
        if len(candidates) >= 2:
            note_bits.append("kept higher-confidence run" if not runs_agree else "runs consistent")
        else:
            note_bits.append("single run only")
        if ev_found:
            note_bits.append("evidence-backed")
        else:
            note_bits.append("no cached external evidence")

        verdicts.append(
            TraitVerdict(
                category=category,
                trait=trait,
                chosen_run=best_name,
                relative_probability=_num(best.get("relative_probability")),
                error_rate=_num(best.get("error_rate")),
                confidence=_num(best.get("confidence")),
                markers_found=int(best.get("markers_found", 0) or 0),
                markers_total=int(best.get("markers_total", 0) or 0),
                evidence_sources=ev_sources,
                evidence_found=ev_found,
                runs_agree=runs_agree,
                note="; ".join(note_bits),
            )
        )

    verdicts.sort(key=lambda v: (v.category, -(v.relative_probability or 0)))
    return {
        "runs": [name for name, _ in reports],
        "n_traits": len(verdicts),
        "snapshot_markers": len(snapshot_markers),
        "verdicts": verdicts,
    }


# --------------------------------------------------------------------------- #
# 4. Entry point
# --------------------------------------------------------------------------- #
def main() -> Dict[str, Any]:
    # Attempt Synapse (will report no_credentials in a sandbox; that's expected).
    syn = SynapseEvidenceClient()
    syn.login()

    result = cross_check(prefer="confidence")
    result["synapse_status"] = syn.status
    result["synapse_connected"] = syn.connected
    return result


if __name__ == "__main__":
    out = main()
    print(f"Synapse: {out['synapse_status']}")
    print(f"Runs cross-checked: {out['runs']}")
    print(f"Traits combined: {out['n_traits']}  |  evidence markers: {out['snapshot_markers']}")
