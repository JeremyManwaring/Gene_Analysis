# Gene_Analysis
Drop your genome.txt file in the folder and run any script.
No setup required: Scripts auto-install pandas if needed, with no VM


Modules include:
**Ethnicity/Ancestry**: Population markers, trait breakdown, interpretation.

**Disease** Risk: Major disease SNPs, polygenic risk scores, clear summaries.

**Longevity**: Aging and protective markers, risk scores, recommendations.

**Fitness**: Athletic trait probabilities and insights.

**Note: These files only look at the existing genes documented by researchers. They do not medically diagnose or give an accurate estimate on genetic variations due to 
possible factors not yet discovered. These datasets also disproportionately draw from European datasets and do not accurately represent every population.**

## Required input

All scripts expect a genotype file named `Genome.txt` in the project
directory. This file should be a tab-separated text file containing four
columns in the following order:

```
rsid    chromosome  position  genotype
```

Lines beginning with `#` are treated as comments and are ignored. The
format mirrors the raw data files provided by personal genomics services
(such as 23andMe).

### Obtaining a genome file

To use your own data, download the raw DNA data from your provider and
rename the file to `Genome.txt` before running any of the analysis
scripts. 23andMe users can get this file by logging in, choosing
**Settings → DNA Relatives → Download Raw Data**, and saving the resulting
text file.

### Example `Genome.txt`

Below is a short sample illustrating the expected format:

```
# rsid    chromosome  position  genotype
rs3094315	1	742429	AA
rs12562034	1	758311	GG
rs3934834	1	995669	CT
```

Place your own `Genome.txt` in the repository root before running any of
the scripts.

_For education only. Not medical advice._ 
