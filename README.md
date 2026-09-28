# Gene Analysis

**Check a raw DNA file against current public genetic evidence, locally.**

One command checks your input, matches measured variants to ClinVar, adds ClinGen gene–disease context, and writes a readable report with a reproducibility record. Your DNA file stays on your computer.

**This is a research annotation tool, not a clinically validated test.** It does not calculate personal disease percentages, predict when a condition will appear, or diagnose a condition. The previous heuristic percentages were removed because they were not validated probabilities.

## Start here

You need **Python 3.11 or newer**. Download or clone this repository, open a terminal in its folder, and install it:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install .
```

On Windows, use `python` instead of `python3`, and activate with `.venv\Scripts\Activate.ps1` in PowerShell.

First, try the entire flow with fictional data. This needs no internet connection after installation:

```sh
gene-analysis demo
```

Then analyze your raw DNA file—no renaming or moving it is needed:

```sh
gene-analysis analyze "/path/to/your/dna.txt"
```

The first real run downloads public reference data (ClinVar is approximately 200 MB compressed and needs additional disk space for its local database). Later runs reuse the cache and refresh it when the source release exceeds the configured age, seven days by default. An update failure is an error, not a successful refresh.

If the file does not say which reference build and strand it uses, confirm those with your provider. For a file confirmed to use GRCh37 and the reference plus strand:

```sh
gene-analysis analyze "/path/to/your/dna.txt" --assembly GRCh37 --strand plus
```

**Do not guess those settings.** Unknown build or strand means clinical matching is withheld and the report explains why. GRCh38 is also supported. Running `python -m Gene_Analysis` works in place of `gene-analysis`.

## What you get

Open the `.md` report in `analysis_reports/`. The same folder contains complete `.json` results and a `.manifest.json` record for reproducing the run.

The pipeline runs in this order:

1. **Input checks:** format, genome build, strand, missing calls, invalid rows and conflicting duplicates.
2. **ClinVar matching:** exact coordinates and observed alleles in the correct genome build, across the usable input records.
3. **ClinGen context:** published gene–disease validity classifications, including inheritance and evidence strength.
4. **Model readiness:** explicitly records which risk predictions cannot be made from the supplied data.
5. **Reports and provenance:** input/software/database hashes, source dates, coverage, evidence links and limitations.

A ClinVar match means a measured allele matches a database record. Review stars describe the supporting review process; they are not a probability. Conflicting and uncertain classifications stay visible. Gene-level evidence does not prove that a particular variant causes disease. No match does not rule out a condition.

## Supported data

Raw tab- or whitespace-separated files, including gzip-compressed files:

```text
# genome build GRCh37
# strand plus
rsid       chromosome  position  genotype
rs12345    1           123456    AG
```

The five-column form `rsid chromosome position allele1 allele2` is also supported. Two allele columns do not establish biological phasing. Supported genotypes are single-base or two-base A/C/G/T calls. Missing calls, insertion/deletion codes and conflicting duplicate markers are excluded and counted.

**Current clinical matching covers single-nucleotide variants.** It does not interpret indels, structural variants, repeats, copy number, mitochondrial heteroplasmy, haplotypes, recessive phase, medication dosing, or complete gene sequencing. The raw reader rejects VCF/BAM/CRAM; use an appropriate sequencing pipeline for those files. This is not comprehensive carrier screening, ancestry inference or pharmacogenomic interpretation.

## Update or reproduce a run

```sh
gene-analysis update-evidence --assembly GRCh37
gene-analysis database-status --assembly GRCh37
```

The defaults store public reference files in `~/.cache/gene-analysis/`. Use `--cache-dir /path/to/cache` to choose another location. GRCh37 and GRCh38 have separate ClinVar databases. ClinGen supplies the gene–disease table. `database-status` prints saved provenance; it does not check the upstream server.

For an exact replay, keep the same input, software revision, ClinVar manifest **and its referenced SQLite database**, and ClinGen snapshot. The report records their locations/hashes. Pin both sources and disable network access:

```sh
gene-analysis analyze "/path/to/your/dna.txt" --offline \
  --clinvar-manifest "/path/to/clinvar-release.json" \
  --clingen-snapshot "/path/to/clingen-release.json" \
  --output-dir replay
```

Include the original `--assembly` / `--strand` overrides if used. Analysis content is reproducible; timestamps, machine metadata and output paths can differ. Pinned old sources remain old and are visibly labeled. Local imports are labeled user-supplied, never live-verified:

```sh
gene-analysis update-evidence --assembly GRCh37 \
  --source-vcf /path/to/clinvar.vcf.gz --clingen-csv /path/to/clingen.csv
```

## Gene combinations, machine learning and future disease risk

A single genotype file cannot train or validate a disease-prediction model. This repository has no longitudinal cohort, age-specific outcome model, external validation population or clinically calibrated absolute-risk model. It therefore **does not rank your future diseases or fabricate percentages**.

For published polygenic scores, use the existing [PGS Catalog Calculator](https://pgsc-calc.readthedocs.io/en/latest/), which performs score download, allele matching, scoring and quality checks. See [the polygenic workflow guide](docs/POLYGENIC_WORKFLOW.md) for the pinned command and required inputs. A polygenic score is not a personal disease probability, and an additive score does not establish gene interactions.

The former handpicked marker panel is optional (`--include-legacy-catalog`). It is unverified historical context, not the current clinical database. Its old weights and percentages are never used for predictions. Historical GWAS/PGS snapshots can be read with `--snapshot-path`; they are not silently refreshed or treated as validated models.

## Privacy and verification

No genotype upload is part of this pipeline. Internet access downloads public reference files only. Reports contain sensitive genetic information: keep them private. Generated results and common raw-data/database filenames are excluded from Git; check your files before sharing. Existing historical commits may still contain old reports; local ignore rules do not remove Git history.

Run the verification suite:

```sh
python -m pip install '.[test]'
python -m pytest -q
```

The active pipeline uses the Python standard library, with no automatic dependency installation at runtime. CI runs offline tests and the demo on Linux, macOS and Windows, across Python 3.11 and 3.14. See [scientific scope and sources](docs/SCIENTIFIC_SCOPE.md) and [the verification record](docs/VERIFICATION.md).

Old scripts and their historical tests are preserved under `legacy/` for audit. Direct execution is disabled because those heuristic outputs were not clinically calibrated. Use the single entry point above.
