# Published polygenic scores and machine learning

## Existing pipeline

Use the published [PGS Catalog Calculator](https://pgsc-calc.readthedocs.io/en/latest/). The upstream latest release checked on 2026-09-27 was [v2.3.0](https://github.com/PGScatalog/pgsc_calc/releases/tag/v2.3.0). Pin it for reproducibility.

It requires Nextflow and a supported execution environment such as Docker. Install them according to the upstream documentation; Gene Analysis does not silently install them. The upstream synthetic verification command is:

```sh
nextflow run pgscatalog/pgsc_calc -r v2.3.0 -profile test,docker
```

This command is provided for the upstream workflow; it has **not** been run or clinically validated as part of this repository's tests.

For real data, follow the [official input preparation guide](https://pgsc-calc.readthedocs.io/en/latest/how-to/prepare.html) and [samplesheet guide](https://pgsc-calc.readthedocs.io/en/latest/how-to/samplesheet.html). Inputs are suitable VCF, PLINK 1 binary, or PLINK 2 datasets—not a renamed consumer text export. Genome-wide scores may need quality-controlled imputation; do not silently upload a genome to an external imputation service.

A samplesheet for PLINK 2 data has these columns:

```csv
sampleset,path_prefix,chrom,format
cohort,/absolute/path/to/quality_checked_genotypes,,pfile
```

After choosing a published score appropriate to the intended trait and validation population:

```sh
nextflow run pgscatalog/pgsc_calc -r v2.3.0 -profile docker \
  --input samplesheet.csv --target_build GRCh37 \
  --pgs_id PGS_ID_YOU_SELECTED --outdir pgs_results
```

Replace the score placeholder with an actual PGS Catalog identifier; do not use an arbitrary example score as a medical recommendation. Preserve matching/QC reports, score files and their hashes, pipeline revision, software versions, source genotypes and source-build provenance. Do not lower overlap thresholds merely to produce a number. Missing calls are not reference alleles, and scores from different diseases are not comparable disease probabilities.

## Why there is no gene-interaction forecast here

The current repository has a person's genotype file and public evidence annotations. It does not contain an approved training cohort with outcomes, onset ages, censoring, family relatedness, ancestry covariates and follow-up. Public variant classifications are not patient outcome labels. Treating variants in one individual as independent patients would invalidate a machine-learning model.

Before a gene-combination model can report disease risk over time, a research project needs:

- An explicit disease definition and prediction horizon; a defensible cohort and data-use permission.
- Appropriate genotype/phenotype quality control, build/allele harmonization, missingness handling and population structure adjustment.
- Person/family/site-separated training and validation with preprocessing and feature selection fitted only to training data.
- A baseline clinical model and additive polygenic model to test whether interactions improve predictions.
- External validation, discrimination and calibration with uncertainty, subgroup performance and checks for overlap/data leakage.
- For time-to-event estimates: event times, censoring and competing risks, plus transportable baseline risk.
- Clinical validation and expert review for the actual intended use, beyond software tests.

No model was trained and no future-disease ranking was produced during this repository revision. These missing prerequisites are reported explicitly instead of substituting unvalidated numeric predictions.

## Evidence

[Lambert, Wingfield et al. (2024)](https://doi.org/10.1038/s41588-024-01937-x) describe the reproducible PGS Catalog calculation workflow. The [eMERGE implementation study](https://doi.org/10.1038/s41591-024-02796-z) describes a separate clinical implementation and validation process, including population calibration; installing a scoring tool does not reproduce that validation.
