"""Conservative reader for 4/5-column raw SNP text; not a VCF/phasing tool."""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional

VALID_BASES = {"A", "C", "G", "T"}
NO_CALLS = {"", "-", "--", "0", "00", "N", "NN", ".", "..", "./.", ".|."}
_ASSEMBLY_ALIASES = {"grch37": "GRCh37", "hg19": "GRCh37", "37": "GRCh37",
                     "grch38": "GRCh38", "hg38": "GRCh38", "38": "GRCh38"}
_MARKER_ID = re.compile(r"(?:rs|i)[0-9]+", re.IGNORECASE)


def _normalize_genotype(genotype: str) -> Optional[str]:
    """Keep haploid calls haploid and normalize unphased diploid A/C/G/T calls."""
    gt = genotype.strip().upper()
    if len(gt) not in (1, 2) or any(base not in VALID_BASES for base in gt):
        return None
    return "".join(sorted(gt))


def _normalize_phased(a1: str, a2: str) -> Optional[str]:
    """Legacy helper name: two raw allele columns do not establish phasing."""
    left, right = a1.strip().upper(), a2.strip().upper()
    if left not in VALID_BASES or right not in VALID_BASES:
        return None
    return "".join(sorted(left + right))


@dataclass
class GenomeLoadStats:
    """Counts describe unique usable markers and excluded input rows."""
    total_lines: int = 0
    parsed_rows: int = 0
    skipped_comment_or_empty: int = 0
    skipped_invalid_rows: int = 0
    skipped_invalid_genotype: int = 0
    inferred_format: str = "unknown"
    skipped_header_rows: int = 0
    skipped_no_calls: int = 0
    skipped_unsupported_genotypes: int = 0
    duplicate_rows: int = 0
    conflicting_rows: int = 0
    excluded_conflicting_markers: int = 0
    excluded_conflicting_loci: int = 0
    haploid_calls: int = 0


class GenomeData:
    """Genotype lookups plus explicit coordinate and strand provenance."""
    def __init__(self, genotype_map: Dict[str, str], records: Optional[List[Dict[str, str]]] = None,
                 stats: Optional[GenomeLoadStats] = None, *, assembly: str = "unknown",
                 assembly_source: str = "unknown", strand: str = "unknown",
                 strand_source: str = "unknown", provider: str = "unknown"):
        self._rsid_index = genotype_map
        self._records = records or []
        self._stats = stats or GenomeLoadStats()
        self.assembly = assembly
        self.assembly_source = assembly_source
        self.strand = strand
        self.strand_source = strand_source
        self.provider = provider

    def get_genotype(self, rsid: str) -> Optional[str]:
        return self._rsid_index.get(rsid)

    def has_snp(self, rsid: str) -> bool:
        return rsid in self._rsid_index

    def get_all_genotypes(self) -> Dict[str, str]:
        return dict(self._rsid_index)

    def get_dataframe(self) -> List[Dict[str, str]]:
        """Legacy row-list helper; copies protect the stored input records."""
        return [dict(record) for record in self._records]

    @property
    def stats(self) -> GenomeLoadStats:
        return self._stats


def _iter_lines(path: Path) -> Iterable[str]:
    # Inspect compressed bytes rather than trusting the filename alone.
    with path.open("rb") as probe:
        compressed = probe.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else Path.open
    with opener(path, "rt", encoding="utf-8-sig", errors="strict") as handle:
        for line in handle:
            yield line.rstrip("\r\n")


def _normalize_assembly(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    normalized = _ASSEMBLY_ALIASES.get(str(value).strip().lower())
    if normalized is None:
        raise ValueError("Supported genome assemblies are GRCh37/hg19 or GRCh38/hg38.")
    return normalized


def _header_assemblies(line: str) -> set[str]:
    values = set()
    for match in re.finditer(r"\b(?:GRCh|hg|NCBI\s*|build\s*)(\d+)\b", line, re.I):
        prefix, number = match.group(0).lower(), match.group(1)
        values.add(_normalize_assembly("hg" + number if prefix.startswith("hg") else number))
    for match in re.finditer(r"\bgenome\s+(?:reference\s+)?version\s+(\d+)(?:\.\d+)?", line, re.I):
        values.add(_normalize_assembly(match.group(1)))
    return values


def _header_strands(line: str) -> set[str]:
    found = set()
    for label, words, sign in (("plus", "plus|positive", r"\+"), ("minus", "minus|negative", "-")):
        if re.search(rf"\b(?:{words})[ -]+(?:\({sign}\)\s*)?strand\b", line, re.I) or re.search(
            rf"\bstrand(?:\s*[:=]\s*|\s+)(?:{words}|{sign})(?:\s|$)", line, re.I
        ):
            found.add(label)
    return found


def _chromosome(value: str, provider: str) -> Optional[str]:
    chrom = value.strip().upper().removeprefix("CHR")
    if chrom in {"M", "MTDNA"}:
        return "MT"
    if provider == "AncestryDNA":
        chrom = {"23": "X", "24": "Y", "25": "PAR", "26": "MT"}.get(chrom, chrom)
    if chrom in {"X", "Y", "MT", "PAR"}:
        return chrom
    if chrom.isdigit() and 1 <= int(chrom) <= 22:
        return str(int(chrom))
    return None


def _resolve_provenance(explicit: Optional[str], headers: set[str], label: str) -> tuple[str, str]:
    if len(headers) > 1:
        raise ValueError(f"Conflicting {label} declarations in the input header.")
    header = next(iter(headers), None)
    if explicit is not None and header is not None and explicit != header:
        raise ValueError(f"Requested {label} {explicit} conflicts with the input header ({header}).")
    if explicit is not None:
        return explicit, "user_provided_and_header" if header else "user_provided"
    return (header, "input_header") if header else ("unknown", "unknown")


def load_genome_data(file_path: str = "Genome.txt", validate: bool = True, *,
                     assembly: Optional[str] = None, strand: Optional[str] = None) -> GenomeData:
    """Read tab/whitespace separated 4/5-column SNP text, optionally gzip.

    Assembly/strand are never inferred from marker positions or genotypes.
    Overrides must agree with explicit headers. No-calls/indels are excluded.
    Conflicting duplicate IDs are quarantined. Haploid X/Y/MT calls are retained
    without manufacturing a second allele. ``validate`` remains for call
    compatibility; empty/unsupported inputs always fail, including when false.
    """
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Genome file not found: {file_path}. Provide a raw genotype text file.")
    if re.search(r"\.vcf(?:\.gz)?$", path.name, re.I):
        raise ValueError("VCF input is not supported. Use a 4/5-column consumer raw SNP text export.")
    explicit_assembly = _normalize_assembly(assembly)
    if strand is not None:
        strand = strand.strip().lower()
        if strand not in {"plus", "minus"}:
            raise ValueError("Strand must be plus or minus; omit it when unknown.")
    stats = GenomeLoadStats()
    by_id: Dict[str, Dict[str, str]] = {}
    locus_calls: dict[tuple[str, str], set[str]] = {}
    conflicts: set[str] = set()
    observed_formats = set()
    header_assemblies: set[str] = set()
    header_strands: set[str] = set()
    provider = "unknown"
    for line in _iter_lines(path):
        stats.total_lines += 1
        stripped = line.strip()
        if stripped.lower().startswith("##fileformat=vcf") or re.match(r"#?CHROM\s+POS\s+ID\s+REF\s+ALT\b", stripped, re.I):
            raise ValueError("VCF input is not supported. Use a 4/5-column consumer raw SNP text export.")
        if not stripped or stripped.startswith("#"):
            stats.skipped_comment_or_empty += 1
            if re.search(r"\bAncestryDNA\b", stripped, re.I):
                provider = "AncestryDNA"
            elif "23andme" in stripped.lower():
                provider = "23andMe"
            header_assemblies.update(_header_assemblies(stripped))
            header_strands.update(_header_strands(stripped))
            continue
        # Keep empty tab fields: a missing allele must not become a haploid call.
        parts = [part.strip() for part in (line.split("\t") if "\t" in line else stripped.split())]
        if [part.lower() for part in parts[:3]] in (["rsid", "chromosome", "position"], ["rsid", "chromosome", "pos"]):
            stats.skipped_header_rows += 1
            continue
        if len(parts) > 5:
            raise ValueError(f"Unsupported {len(parts)}-column input on line {stats.total_lines}; expected 4 or 5 columns.")
        if len(parts) not in {4, 5}:
            stats.skipped_invalid_rows += 1
            continue
        marker, chrom = parts[0].lower(), _chromosome(parts[1], provider)
        if _MARKER_ID.fullmatch(marker) is None or chrom is None or not parts[2].isdigit() or int(parts[2]) < 1:
            stats.skipped_invalid_rows += 1
            continue
        observed_formats.add(len(parts))
        no_call = any(allele.upper() in NO_CALLS for allele in parts[3:])
        genotype = _normalize_genotype(parts[3]) if len(parts) == 4 else _normalize_phased(parts[3], parts[4])
        if genotype is None or (len(genotype) == 1 and chrom not in {"X", "Y", "MT"}):
            stats.skipped_invalid_genotype += 1
            if no_call:
                stats.skipped_no_calls += 1
            else:
                stats.skipped_unsupported_genotypes += 1
            continue
        record = {"rsid": marker, "chromosome": chrom, "position": str(int(parts[2])), "genotype": genotype}
        locus_calls.setdefault((chrom, record["position"]), set()).add(genotype)
        if marker in conflicts:
            stats.duplicate_rows += 1
            stats.conflicting_rows += 1
            continue
        previous = by_id.get(marker)
        if previous is not None:
            stats.duplicate_rows += 1
            if previous != record:
                stats.conflicting_rows += 1
                conflicts.add(marker)
                del by_id[marker]
            continue
        by_id[marker] = record
    resolved_assembly, assembly_source = _resolve_provenance(explicit_assembly, header_assemblies, "genome assembly")
    resolved_strand, strand_source = _resolve_provenance(strand, header_strands, "strand")
    conflicting_loci = {locus for locus, calls in locus_calls.items() if len(calls) > 1}
    # Separate probes/identifiers cannot provide contradictory genotypes at one
    # coordinate. Quarantine the entire locus rather than selecting a probe.
    for marker, row in list(by_id.items()):
        if (row["chromosome"], row["position"]) in conflicting_loci:
            conflicts.add(marker)
            del by_id[marker]
    records = list(by_id.values())
    stats.parsed_rows = len(records)
    stats.excluded_conflicting_markers = len(conflicts)
    stats.excluded_conflicting_loci = len(conflicting_loci)
    stats.haploid_calls = sum(len(row["genotype"]) == 1 for row in records)
    stats.inferred_format = {frozenset({4}): "4col_genotype", frozenset({5}): "5col_allele1_allele2",
                             frozenset({4, 5}): "mixed_4col_5col"}.get(frozenset(observed_formats), "unknown")
    if not records:
        raise ValueError("No valid, nonconflicting genotype rows were parsed. Expected 4/5-column raw SNP text with A/C/G/T calls.")
    return GenomeData(genotype_map={marker: row["genotype"] for marker, row in by_id.items()},
                      records=records, stats=stats, assembly=resolved_assembly, assembly_source=assembly_source,
                      strand=resolved_strand, strand_source=strand_source, provider=provider)
