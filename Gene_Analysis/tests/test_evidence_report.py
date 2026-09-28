"""Evidence reports must distinguish unknown results from measured observations."""

from datetime import UTC, datetime, timedelta

from Gene_Analysis.analysis.catalog import MarkerDefinition
from Gene_Analysis.analysis.pipeline import build_analysis_report
from Gene_Analysis.utils.genome_loader import GenomeData


def marker(rsid="rs1", trait="Type 2 diabetes"):
    return MarkerDefinition("disease", trait, rsid, "T", "TEST", "Historical catalog label")


def snapshot():
    return {
        "metadata": {"created_at": datetime.now(UTC).isoformat()},
        "markers": {
            "rs1": {
                "best_last_updated": datetime.now(UTC).isoformat(),
                "gwas_associations": [{
                    "association_id": "test-association", "reported_trait": ["Type 2 diabetes"],
                    "snp_allele": [{"rs_id": "rs1", "effect_allele": "T"}],
                    "or_value": 1.4, "p_value": 1e-10, "ci_lower": 1.2, "ci_upper": 1.6,
                }],
            }
        },
        "trait_models": {"Type 2 diabetes": {"performance": {"metrics": {
            "class_acc": [{"name_short": "AUROC", "estimate": 0.9}]
        }}}},
    }


def test_no_prediction_in_rows_or_aggregates_even_with_study_evidence() -> None:
    report = build_analysis_report(GenomeData({"rs1": "TT"}), snapshot(), [marker()])
    for section in ("rows", "trait_summaries", "category_summaries"):
        for row in report[section]:
            for key in ("relative_probability", "error_rate", "confidence", "heritability_proxy"):
                assert row[key] is None
    row = report["rows"][0]
    assert row["study_odds_ratio"] == 1.4
    assert row["study_beta"] is None
    assert row["effect_size"] is None
    assert row["raw_effect_allele_count"] == 2
    assert row["risk_allele_count"] is None
    assert row["allele_harmonization"] == "not_performed"
    assert report["model_status"]["gene_combination_machine_learning"] == "not_trained_or_validated"
    assert report["model_status"]["longitudinal_disease_probabilities"] == "not_estimated"
    assert report["source_evidence"]["markers"]["rs1"]["gwas_associations"]


def test_missing_genotype_does_not_look_protective_or_homozygous_reference() -> None:
    report = build_analysis_report(GenomeData({}), snapshot(), [marker()])
    row = report["rows"][0]
    assert row["genotype_status"] == "missing"
    assert row["raw_effect_allele_count"] is None
    assert row["risk_allele_count"] is None
    assert report["metadata"]["markers_found"] == 0


def test_catalog_allele_never_substitutes_for_unverified_source_allele() -> None:
    evidence = snapshot()
    evidence["markers"]["rs1"]["gwas_associations"][0]["snp_allele"] = []
    row = build_analysis_report(GenomeData({"rs1": "TT"}), evidence, [marker()])["rows"][0]
    assert row["catalog_allele"] == "T"
    assert row["effect_allele"] is None
    assert row["study_odds_ratio"] is None
    assert row["association_status"] == "effect_allele_unverified_or_mismatched"


def test_source_edit_date_does_not_imply_recent_retrieval() -> None:
    report = build_analysis_report(GenomeData({}), snapshot(), [marker()])
    assert report["rows"][0]["evidence_freshness"] == "unknown"
    evidence = snapshot()
    evidence["markers"]["rs1"]["retrieved_at"] = (datetime.now(UTC) - timedelta(days=31)).isoformat()
    report = build_analysis_report(GenomeData({}), evidence, [marker()])
    assert report["rows"][0]["evidence_freshness"] == "stale"


def test_source_errors_do_not_count_as_annotations() -> None:
    evidence = {"markers": {"rs1": {"clinvar": {"error": "network failure"}, "gnomad": {"skipped": True}}}}
    row = build_analysis_report(GenomeData({"rs1": "TT"}), evidence, [marker()])["rows"][0]
    assert row["source"] == ""


def test_explicitly_empty_catalog_stays_empty() -> None:
    report = build_analysis_report(GenomeData({"rs1": "TT"}), {}, [])
    assert report["rows"] == []
    assert report["trait_summaries"] == []
    assert report["metadata"]["total_markers"] == 0


def test_unique_marker_coverage_does_not_double_count_trait_labels() -> None:
    report = build_analysis_report(GenomeData({"rs1": "TT"}), {}, [marker(), marker(trait="Other")])
    assert report["metadata"]["total_markers"] == 2
    assert report["metadata"]["unique_catalog_markers"] == 1
    assert report["metadata"]["unique_catalog_markers_found"] == 1


def clinical_report():
    report = build_analysis_report(GenomeData({"rs1": "CT"}), {}, [])
    report["clinical"] = {
        "status": "complete", "scope": "Only exact supported SNVs are matched.",
        "database": {"source_release": "20260101", "source_sha256": "abc123", "stale": True},
        "summary": {"input_records": 1, "matched_variants": 1, "professional_review_candidates": 1},
        "findings": [{
            "gene": "TEST:1", "chromosome": "1", "position": 100,
            "reference": "C", "alternate": "T", "genotype": "CT",
            "germline_classification": "Pathogenic", "review_status": "reviewed_by_expert_panel",
            "review_stars": 3, "needs_professional_review": True,
            "conditions": "Example | condition", "clinvar_variation_id": "123",
            "absolute_risk": None, "time_dependent_risk": None,
        }],
    }
    return report


def test_human_report_explains_source_classifications_and_unknown_predictions(tmp_path) -> None:
    from Gene_Analysis.analysis.reporting import write_markdown_report
    report = clinical_report()
    report["metadata"]["synthetic_demo"] = True
    report["metadata"]["pipeline_steps"] = [{"step": "input_validation", "status": "completed"}]
    output = tmp_path / "report.md"
    write_markdown_report(report, str(output))
    text = output.read_text()
    assert "SYNTHETIC DEMO" in text
    assert "Pathogenic" in text
    assert "3/4 stars" in text
    assert "Synthetic fixture; no live record" in text
    assert "https://www.ncbi.nlm.nih.gov/clinvar/variation/123/" not in text
    assert "20260101" in text
    assert "abc123" in text
    assert "age threshold: **yes**" in text
    assert "not a clinically validated" in text
    assert "not a personal confidence score" in text
    assert "Example \\| condition" in text
    assert "| Relative Probability |" not in text
    assert "| Error Rate |" not in text
    assert "Heritability Proxy" not in text
    assert "input validation" in text


def test_missing_clinical_matching_is_visibly_withheld(tmp_path, capsys) -> None:
    from Gene_Analysis.analysis.reporting import print_console_summary, write_markdown_report
    report = build_analysis_report(GenomeData({}), {}, [])
    report["clinical"] = {"status": "withheld", "reason": "Unknown genome build."}
    output = tmp_path / "withheld.md"
    write_markdown_report(report, str(output))
    text = output.read_text()
    assert "Clinical matching was not completed" in text
    assert "Unknown genome build." in text
    assert "negative clinical test" in text
    print_console_summary(report)
    assert "Clinical matching withheld: Unknown genome build." in capsys.readouterr().out


def test_markdown_truncation_is_explicit_and_full_json_keeps_all_matches(tmp_path) -> None:
    import json
    from Gene_Analysis.analysis.reporting import write_json_report, write_markdown_report
    report = clinical_report()
    template = report["clinical"]["findings"][0]
    report["clinical"]["findings"] = [dict(template, clinvar_variation_id=str(i)) for i in range(105)]
    output = tmp_path / "all.md"
    write_markdown_report(report, str(output))
    assert "Showing 100 of 105" in output.read_text()
    json_output = tmp_path / "all.json"
    write_json_report(report, str(json_output))
    assert len(json.loads(json_output.read_text())["clinical"]["findings"]) == 105


def test_catalog_console_order_is_alphabetical_not_legacy_confidence(capsys) -> None:
    from Gene_Analysis.analysis.reporting import print_console_summary
    report = build_analysis_report(GenomeData({}), {}, [])
    report["trait_summaries"] = [
        {"trait": "Z condition", "confidence": 99, "relative_probability": 99},
        {"trait": "A condition", "confidence": 1, "relative_probability": 1},
    ]
    print_console_summary(report)
    output = capsys.readouterr().out
    assert output.index("A condition") < output.index("Z condition")
    assert "99%" not in output
    assert "1%" not in output


def test_clingen_database_scope_and_matched_context_are_distinct(tmp_path) -> None:
    from Gene_Analysis.analysis.reporting import write_markdown_report
    report = clinical_report()
    report["clinical"]["findings"][0]["clingen_gene_validity"] = [{
        "gene_symbol": "TEST", "classification": "Definitive",
        "disease_identifier_in_clinvar_aggregate": False,
    }]
    report["gene_validity"] = {
        "status": "annotated", "snapshot_sha256": "test-hash", "cache_stale": False,
        "metadata": {
            "source_created_at": "2026-01-01", "fetched_at": "2026-01-02",
            "summary": {"curations": 10, "genes": 8, "strong_or_definitive_curations": 6,
                        "genes_with_strong_or_definitive_curation": 5,
                        "classification_counts": {"Definitive": 6, "Refuted": 4}},
        },
    }
    output = tmp_path / "context.md"
    write_markdown_report(report, str(output))
    text = output.read_text()
    assert "**10** gene–disease curations" in text
    assert "**8** genes" in text
    assert "**1 / 1**" in text
    assert "**0** have an exact disease identifier" in text
    assert "| Refuted | 4 |" in text
    assert "Source export created: `2026-01-01`" in text
    assert "Snapshot SHA-256: `test-hash`" in text
