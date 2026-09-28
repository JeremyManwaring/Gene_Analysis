# Verification record

Completed 2026-09-27 (America/Los_Angeles). This verifies software behavior and source integrity. **It is not clinical validation.**

## Checks performed

- **107 automated tests passed** in the working checkout and again in a clean export of the candidate Git tree.
- Installed a normal wheel into a new environment, then ran the offline demo from outside the repository with isolated Python imports. The installed package completed the entire annotation/report flow without depending on the original checkout.
- The export contains ordinary source files, with **zero Git submodule/gitlink entries**. This fixes the previous nested checkout that a fresh clone could not recreate.
- The core runtime has no third-party Python dependencies. Python 3.14.7/macOS were exercised locally; Python 3.11/3.14 on Linux/macOS/Windows are configured in CI but those remote jobs have **not** been run in this session.
- Input tests cover raw four/five-column text, gzip, explicit build/strand, unknown or contradictory provenance, haploid calls, malformed rows, no-calls, duplicated identifiers and conflicting calls at one coordinate.
- Evidence tests cover exact alternate-allele matching, reference-only calls, wrong build/strand, unsupported variants, review stars, uncertain/conflicting classifications and multiple submitters with **no conflicts**.
- Integrity checks cover modified snapshots/databases, failed refresh preservation, incomplete/old evidence, source-based freshness, a completed interrupted download, and the ISO-formatted file date used by the actual NCBI release.
- The network-disabled demo produced reproducible matching content. Missing pinned files and invalid inputs fail clearly; outputs cannot overwrite the input file.
- A full live-database smoke test matched **30 of 30 synthetic input loci** against the official ClinVar index, attached ClinGen context and wrote JSON/Markdown/manifests. These are constructed test calls, not personal results.
- Personal disease probabilities, confidence/error percentages, heritability substitutes, gene-interaction forecasts and disease-onset rankings are not produced. Tests enforce this even when old source statistics are supplied.

## Live reference snapshot

| Item | Verified value |
|---|---|
| ClinVar assembly | GRCh37 |
| ClinVar source file date | 2026-09-23 |
| NCBI Last-Modified | Thu, 24 Sep 2026 09:06:57 GMT |
| Full VCF records indexed | 4,555,389 |
| Single-nucleotide records | 4,234,744 |
| Compressed source size | 198,082,810 bytes |
| Upstream published MD5 matched | Yes |
| SQLite integrity check | Passed |
| ClinGen source date | 2026-09-27 |
| ClinGen gene–disease curations | 3,679 |
| Distinct HGNC genes in those curations | 3,038 |
| Strong/Definitive curations | 2,379 |

ClinVar source SHA-256:
`2990cf7124942759c417e22fe31fb8eb668b2797cd93fb8f924b9b9b2a1ad8a3`

ClinVar database SHA-256:
`53ec83e3eb61a62f52b535695572f40802042bef8b306d8a6c09c392e6fadba8`

ClinGen source SHA-256:
`ecabebf910705adc5063fc7f74e6aaf27cf184c1a65682c15c481876884ae6d6`

Sources: [official ClinVar GRCh37 VCF](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz), [official ClinGen validity export](https://search.clinicalgenome.org/kb/gene-validity/download).

The live GRCh37 database was downloaded and verified. GRCh38 uses the same tested importer and has synthetic build-mismatch coverage; the full GRCh38 download was **not** performed during this audit. The default seven-day refresh window is a cache policy, not a claim that every upstream submission is scientifically current.

## Scope that remains unsupported

There is no clinically validated personal-risk or time-to-event model, no outcome-labeled longitudinal training cohort, and no performed machine-learning analysis. The existing PGS Catalog Calculator is documented as the appropriate external score workflow; it was not executed here. This code does not cover all genetic variant classes or every clinically relevant gene.

Old heuristic reports were preserved locally and removed from the current tracked files. Existing Git history was not rewritten. Earlier commits can still contain those reports. The nested repository metadata and pre-migration index are backed up under `.git/local-migration-backup/`.

This record documents the pre-publication audit. Publication status is recorded in Git history.
