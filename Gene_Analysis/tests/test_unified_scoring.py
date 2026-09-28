"""Regression tests against uncalibrated personal risk and unsafe allele matches."""

from __future__ import annotations

import pytest

from Gene_Analysis.analysis.scoring import (
    _safe_float,
    allele_count,
    ancestry_bias_penalty,
    association_effect_allele,
    confidence_from_error_and_bias,
    pgs_trait_model_metrics,
    pick_best_association_for_trait,
    relative_probability_from_effect,
    snp_error_rate,
)


@pytest.mark.parametrize("count", [None, 0, 1, 2])
def test_effects_cannot_be_converted_to_probabilities(count) -> None:
    assert relative_probability_from_effect(count, effect_or=1.6) is None
    assert relative_probability_from_effect(count, beta=1000) is None
    assert relative_probability_from_effect(count, weight=0.8) is None


def test_precision_and_ancestry_are_not_prediction_error_or_confidence() -> None:
    assert snp_error_rate(1.18, 1.22, 1.2, None) is None
    penalty, note = ancestry_bias_penalty({"eval": {"dist": {"EUR": 95, "AFR": 5}}})
    assert penalty is None
    assert "not established" in note
    assert confidence_from_error_and_bias(0.1, 0.1) is None


def test_pgs_metadata_is_preserved_without_personal_score() -> None:
    performance = {
        "metrics": {
            "effect_sizes": [{"estimate": "1.3", "ci_lower": "1.1", "ci_upper": "1.5"}],
            "class_acc": [{"name_short": "AUROC", "estimate": "0.72"}],
            "othermetrics": [{"name_short": "R2", "estimate": "0.18"}],
        }
    }
    result = pgs_trait_model_metrics({"performance": performance})
    assert result["status"] == "metadata_only"
    assert result["published_performance"] == performance
    for key in ("personal_score", "prs_probability", "prs_error_rate", "heritability_proxy", "prs_metric"):
        assert result[key] is None
    assert pgs_trait_model_metrics({"error": "timeout"})["status"] == "unavailable"


@pytest.mark.parametrize("genotype", [None, "", "--", "NN", "A-", "ID", "AAA", "nan"])
def test_missing_or_unsupported_genotype_is_not_reference(genotype) -> None:
    assert allele_count(genotype, "A") is None


def test_literal_counts_preserve_haploidy_and_do_not_complement() -> None:
    assert allele_count("AG", "A") == 1
    assert allele_count("G", "G") == 1
    assert allele_count("AA", "G") == 0
    assert allele_count("C/T", "C") == 1
    assert allele_count("A|A", "A") == 2
    assert allele_count("AA", "AT") is None
    assert allele_count("AA", None) is None


def _association(trait="Type 2 diabetes", allele="T", rsid="rs1", p=1e-9):
    return {
        "reported_trait": [trait], "p_value": p,
        "snp_allele": [{"rs_id": rsid, "effect_allele": allele}],
    }


def test_no_fallback_to_unrelated_or_partially_matching_trait() -> None:
    associations = [_association("Type 1 diabetes", p=1e-100), _association("Diabetes medication use")]
    assert pick_best_association_for_trait(associations, "Type 2 diabetes", "rs1", "T") is None
    assert pick_best_association_for_trait(associations, "Diabetes", "rs1", "T") is None


def test_exact_trait_and_effect_allele_select_correct_study_record() -> None:
    good = _association()
    records = [
        _association("Type 1 diabetes", p=1e-100),
        _association(allele="C", p=1e-50),
        _association(rsid="rs2", p=1e-30), good,
    ]
    assert pick_best_association_for_trait(records, "type-2 diabetes", "rs1", "T") is good


@pytest.mark.parametrize("allele", ["?", "AT", "T,C", "-", ""])
def test_uncertain_or_multiallelic_effect_is_not_assigned(allele) -> None:
    assert association_effect_allele(_association(allele=allele), "rs1") is None


def test_haplotype_effect_is_not_misattributed_to_one_snp() -> None:
    assoc = _association()
    assoc["snp_allele"].append({"rs_id": "rs2", "effect_allele": "A"})
    assert association_effect_allele(assoc, "rs1") is None
    assert pick_best_association_for_trait([assoc], "Type 2 diabetes", "rs1", "T") is None


@pytest.mark.parametrize("value", ["OR 1.3", "NaN", float("inf"), True, "95% CI 1.2-1.4"])
def test_arbitrary_text_and_nonfinite_values_are_not_estimates(value) -> None:
    assert _safe_float(value) is None
