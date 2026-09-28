# Scientific scope and source interpretation

Reviewed against the primary sources below on 2026-09-27. Database freshness is a point-in-time property, not a guarantee that every record reflects every recent publication.

| Source | What it establishes | What it does not establish |
|---|---|---|
| ClinVar full assembly-specific VCF | Submitted variant classifications, conditions and review status for mapped variants | A personal diagnosis, penetrance estimate, complete assay, or independent NIH verification of every claim |
| ClinGen gene–disease validity CSV | Expert-curated strength of a specific gene–disease relationship and inheritance context | That every variant in that gene is pathogenic, or the person's probability of disease |
| Historical GWAS annotations, optional | Population-level statistical associations with study-specific alleles, traits and effects | Personal absolute probability, causality, valid interaction effects or prediction across all ancestries |
| PGS Catalog / pgsc_calc, external | Published scoring definitions and reproducible additive score calculation with QC | Clinical suitability for every user, gene interactions, or a time-dependent probability without calibration |

## Matching and evidence policy

Clinical matching requires an explicitly known assembly, plus-strand calls, and exact chromosome/position/observed alternate allele. The full ClinVar VCF is indexed, but this implementation interprets matching single-nucleotide alleles only. Unmeasured variants are unknown, not negative. Ambiguous or unsupported inputs are excluded and counted. Two allele columns are not evidence of chromosome phase.

All matched classifications remain available in JSON. Expert-panel / practice-guideline review is distinguished from multiple or single submissions. Pathogenic/likely-pathogenic classifications with adequate review can be marked for professional review; they are not diagnoses or confirmed laboratory results. Aggregate variant classifications and condition lists must not be read as a one-to-one disease assertion without reviewing the linked source.

ClinGen classifications include Strong/Definitive alongside Moderate/Limited/Disputed/Refuted and other categories. A gene name match adds context only: it cannot upgrade variant pathogenicity or imply an exact disease match. The public table covers what ClinGen has curated, not every possible clinically relevant gene.

Downloads preserve the source date, transport checksum when supplied, SHA-256, and local database integrity. A failed refresh must not make old or incomplete evidence look current. Immutable pinned files can be used offline, with staleness shown. Checksums demonstrate reproducibility/integrity, not clinical truth.

## Removed calculations

The former pipeline used `OR / (1 + OR)` as a personal-looking percentage, sometimes substituted `exp(beta)`, blended population PGS performance into individual scores, invented error/confidence percentages, and labeled predictive statistics as heritability proxies. A separate legacy module scaled arbitrary baseline prevalence by SNP multipliers and ranked diseases. These methods were not clinically calibrated; none is used by the active pipeline. Historical reports are not retroactively corrected and should not be interpreted as validated predictions.

Exact study odds ratios and beta estimates are different quantities. Their confidence intervals describe the fitted population estimate, not the error probability of this person's result. AUC is not individual accuracy, and prediction R-squared is not interchangeable with heritability. No percentage substitute is generated when required evidence is absent.

## Primary references

- [ClinVar review status and stars](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).
- [ClinVar FTP products and VCF coverage](https://www.ncbi.nlm.nih.gov/clinvar/docs/ftp_primer/).
- [ClinVar data validation](https://www.ncbi.nlm.nih.gov/clinvar/docs/validation/).
- [ClinGen gene–disease validity](https://www.clinicalgenome.org/curation-activities/gene-disease-validity/).
- [ClinGen official bulk files](https://search.clinicalgenome.org/kb/downloads), including the [live gene-validity CSV](https://search.clinicalgenome.org/kb/gene-validity/download).
- [PGS Catalog Calculator documentation](https://pgsc-calc.readthedocs.io/en/latest/) and [publication](https://doi.org/10.1038/s41588-024-01937-x).
- [eMERGE score selection and clinical implementation study](https://doi.org/10.1038/s41591-024-02796-z).

Plugin discovery found the installed Consensus research connector, which was used to find and retrieve the two peer-reviewed pipeline/validation papers. Official source documentation and public downloads supply the database definitions and provenance. No plugin independently certifies this software as clinically validated.
