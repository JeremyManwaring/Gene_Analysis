"""Marker catalog utilities used by evidence update and analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

from ..config import (
    ABO_SNPS,
    ANCESTRY_CATEGORIES,
    CARRIER_SCREENING_SNPS,
    COGNITIVE_SNPS,
    DISEASE_SNPS,
    FITNESS_SNPS,
    LONGEVITY_MARKERS,
    NUTRITION_SNPS,
    PHARMACOGENOMICS_SNPS,
    RHD_SNPS,
    SLEEP_SNPS,
)


@dataclass(frozen=True)
class MarkerDefinition:
    category: str
    trait: str
    rsid: str
    risk_allele: Optional[str]
    gene: Optional[str]
    description: str
    weight: float = 1.0


def _normalize_risk_allele(value: Optional[str]) -> Optional[str]:
    if not value:
        return None
    allele = value.strip().upper()
    if len(allele) == 1 and allele in {"A", "C", "G", "T"}:
        return allele
    return None


def _marker(
    category: str,
    trait: str,
    rsid: str,
    risk_allele: Optional[str],
    description: str,
    gene: Optional[str] = None,
    weight: float = 1.0,
) -> MarkerDefinition:
    return MarkerDefinition(
        category=category,
        trait=trait,
        rsid=rsid,
        risk_allele=_normalize_risk_allele(risk_allele),
        gene=gene,
        description=description,
        weight=weight,
    )


def _gene_from_description(description: str) -> Optional[str]:
    if " - " in description:
        return description.split(" - ", 1)[0].strip()
    return None


def build_marker_catalog() -> List[MarkerDefinition]:
    """Build a unified marker catalog across all feature families."""
    markers: List[MarkerDefinition] = []

    for rsid, data in DISEASE_SNPS.items():
        risk = data.get("risk", [None])[0]
        markers.append(
            _marker(
                category="disease",
                trait=data.get("disease", "Disease Risk"),
                rsid=rsid,
                risk_allele=risk,
                gene=data.get("gene"),
                description=data.get("description", ""),
                weight=1.0,
            )
        )

    for rsid, data in LONGEVITY_MARKERS.items():
        risk = data.get("risk", [None])[0]
        markers.append(
            _marker(
                category="longevity",
                trait="Longevity",
                rsid=rsid,
                risk_allele=risk,
                gene=_gene_from_description(data.get("description", "")),
                description=data.get("description", ""),
                weight=float(data.get("weight", 1.0)),
            )
        )

    for trait, snps in FITNESS_SNPS.items():
        for rsid, data in snps.items():
            risk = data.get("risk", [None])[0]
            markers.append(
                _marker(
                    category="fitness",
                    trait=trait,
                    rsid=rsid,
                    risk_allele=risk,
                    gene=_gene_from_description(data.get("description", "")),
                    description=data.get("description", ""),
                    weight=float(data.get("weight", 1.0)),
                )
            )

    for ancestry_category, snps in ANCESTRY_CATEGORIES.items():
        for rsid, description in snps.items():
            markers.append(
                _marker(
                    category="ancestry",
                    trait=ancestry_category,
                    rsid=rsid,
                    risk_allele=None,
                    gene=_gene_from_description(description),
                    description=description,
                    weight=1.0,
                )
            )

    for rsid, data in PHARMACOGENOMICS_SNPS.items():
        markers.append(
            _marker(
                category="pharmacogenomics",
                trait=", ".join(data.get("drugs", [])) or "Drug Metabolism",
                rsid=rsid,
                risk_allele=data.get("risk_allele"),
                gene=data.get("gene"),
                description=data.get("description", ""),
                weight=1.0,
            )
        )

    for rsid, data in NUTRITION_SNPS.items():
        markers.append(
            _marker(
                category="nutrition",
                trait=data.get("trait", "Nutrition"),
                rsid=rsid,
                risk_allele=data.get("risk_allele"),
                gene=data.get("gene"),
                description=data.get("description", ""),
                weight=1.0,
            )
        )

    for rsid, data in SLEEP_SNPS.items():
        markers.append(
            _marker(
                category="sleep",
                trait=data.get("trait", "Sleep"),
                rsid=rsid,
                risk_allele=data.get("risk_allele"),
                gene=data.get("gene"),
                description=data.get("description", ""),
                weight=1.0,
            )
        )

    for rsid, data in ABO_SNPS.items():
        markers.append(
            _marker(
                category="blood_type",
                trait="ABO",
                rsid=rsid,
                risk_allele=None,
                gene="ABO",
                description=data.get("description", ""),
                weight=1.0,
            )
        )

    for rsid, data in RHD_SNPS.items():
        markers.append(
            _marker(
                category="blood_type",
                trait="Rh",
                rsid=rsid,
                risk_allele=None,
                gene="RHD",
                description=data.get("description", ""),
                weight=1.0,
            )
        )

    for rsid, data in CARRIER_SCREENING_SNPS.items():
        markers.append(
            _marker(
                category="carrier_screening",
                trait=data.get("condition", "Carrier Screening"),
                rsid=rsid,
                risk_allele=data.get("risk_allele"),
                gene=data.get("gene"),
                description=data.get("description", ""),
                weight=1.0,
            )
        )

    for domain, snps in COGNITIVE_SNPS.items():
        for rsid, data in snps.items():
            markers.append(
                _marker(
                    category="cognitive",
                    trait=domain,
                    rsid=rsid,
                    risk_allele=data.get("risk_allele"),
                    gene=data.get("gene"),
                    description=data.get("description", ""),
                    weight=1.0,
                )
            )

    return markers


def marker_groups(markers: Iterable[MarkerDefinition]) -> Dict[str, List[MarkerDefinition]]:
    grouped: Dict[str, List[MarkerDefinition]] = {}
    for marker in markers:
        grouped.setdefault(marker.category, []).append(marker)
    return grouped


def marker_rsids(markers: Iterable[MarkerDefinition]) -> List[str]:
    return sorted({marker.rsid for marker in markers})
