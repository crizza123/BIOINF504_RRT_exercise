## The Reproducibility Audit: Find the Flaws, Fix the Pipeline

**Estimated time:** ~60 minutes  
**Format:** Individual  
**Deliverable:** Completed tasks submitted as a pull request or a single PDF or `.md` file

#### Created with the use of claude.ai

---

### Learning Objectives

By the end of this exercise you should be able to:

- Identify p-hacking and underpowered study design
- Spot missing randomization and confounded replicates
- Evaluate code for versioning and reproducibility gaps
- Apply FAIR principles to data and metadata
- Recognize sampling bias and demographic confounders
- Propose concrete, actionable remediation steps
- Prioritize which flaws matter most in a real-world setting

---

### Scenario

A collaborator shares a GitHub repository for a published RNA-seq differential expression study comparing gene expression in tumor vs. normal tissue across two cohorts. They are excited — several genes hit *p* < 0.05 and they want to submit a follow-up grant. Before you endorse the work, you agree to do a quick reproducibility audit.

**Your job: find everything that could go wrong, and recommend fixes.**

---

### Background: Key Terms

Before you begin, make sure you are comfortable with the following terms used throughout this assignment.

- **HARKing** (Hypothesizing After Results are Known): reporting a finding as if it were an a priori hypothesis when it was actually identified after inspecting the data. This inflates the apparent significance of results.
- **Multiple testing correction**: when testing many hypotheses simultaneously (e.g., 18,000 genes), the probability of at least one false positive at α = 0.05 approaches certainty. Procedures such as Benjamini–Hochberg FDR control the expected proportion of false discoveries.
- **Negative binomial model**: RNA-seq counts are discrete, non-negative integers. They are overdispersed relative to a Poisson distribution (variance > mean), which makes the negative binomial distribution a better fit. Tools like DESeq2 (R) and pydeseq2 (Python) use this model.
- **CV (coefficient of variation)**: CV = σ/μ, the ratio of standard deviation to mean. A CV of 0.4 means the standard deviation is 40% of the mean — a typical level of biological variability in RNA-seq data.
- **FAIR principles**: a framework for data sharing — data should be Findable, Accessible, Interoperable, and Reusable.
- **Batch effect**: systematic technical variation introduced by differences in sample processing across time, lab, reagent lot, or sequencing run, which can masquerade as biological signal.

---

### Step 1 — Read the Methods Snapshot *(~7 min)*

#### Provided methods excerpt

> "We analyzed RNA-seq data from Cohort A (n = 6 tumor, n = 6 normal) and Cohort B (n = 4 tumor, n = 4 normal). All patients were male, aged 55–70, recruited at a single institution between 2018–2022. Library preparation and sequencing were performed in two runs: samples collected in 2018 were processed together, and samples collected in 2022 were processed together. Three tumor samples from Patient_01 were re-sequenced at higher depth and included as separate entries. We tested 18,000 genes for differential expression using a t-test and reported all genes with p < 0.05 without multiple testing correction. We selected the top 12 candidate genes based on fold-change ranking after seeing the data. The analysis was run in [Python/R]; the version is not recorded and no environment file was committed. Raw FASTQ files are available upon reasonable request. The pipeline script is included in the repo but requires manual setup with no version guidance. All samples are labeled 'Patient_01' through 'Patient_20' with no additional metadata file. Both authors are from the same lab and share the same clinical specialty (oncology)."

#### Task 1A — Annotation checklist

For each item below, mark whether a problem is present and write one sentence describing it. If a problem is not present or cannot be determined from the text, say so and explain.

| # | Category | Problem present? | Description |
|---|---|---|---|
| 1 | Sample size / statistical power | Yes | Cohort A and B are below the 20-25 group requirement for differential gene detection. |
| 2 | Multiple testing correction | Yes | 18000 genes are tested with no correction.|
| 3 | P-hacking / HARKing | Yes | The filter threshold was choosen after looking at the results.|
| 4 | Appropriate statistical test for count data | Yes | A two-sample t-test is not appropriate for RNA-seq because RNA-seq counts are discreet and overdispersed.|
| 5 | Randomization and batch effects |Yes | The year of sample prep is a confounding variable because Cohort B is all done in 2022. |
| 6 | Biological vs. technical replicates |Yes | Patient 1 was sequenced three times but each library is counted as a biological replicate. |
| 7 | Version pinning and dependency management |Yes | There is no version control, env.yml, or requirements.txt. |
| 8 | Code portability / containerization | Yes | There is no contanerization through docker or singularity. |
| 9 | FAIR — Findability and Accessibility of raw data | Yes | No access to the raw data or upload to GEO / other repositories for public to use. |
| 10 | FAIR — Interoperability and metadata completeness | Yes | There is no metadata file to use so we do not know how it relates to the samples. |
| 11 | FAIR — Reusability, licensing, and provenance | Yes | No license file in the repository. |
| 12 | Single-institution / demographic sampling bias | Yes| All patients are male aged 55-70, no females included biases the results. |
| 13 | Confirmation bias from homogeneous team | Yes | Both authors are from the same lab and specialty. |
| 14 | Domain blind spots from homogeneous team | Yes | There are no domain experts in bioinformatics or statistics. |

> **Hint:** Not every category has a single flaw — some have multiple layered problems. The goal is systematic thinking, not a perfect list.

---

### Step 2 — Inspect the Code *(~12 min)*

Choose **Python** or **R** below. Review the script for software engineering and reproducibility issues.

---

#### Option A — Python

```python
# analysis.py
# (no author, no date, no Python version, no docstring)

import pandas as pd
import numpy as np
from scipy import stats
import os

os.chdir("/Users/mitreacristina/Desktop/rnaseq_project")

counts   = pd.read_csv("counts_final_FINAL_v3.csv", index_col=0)
metadata = pd.read_csv("meta.csv")

# filter low counts — threshold chosen after looking at results
counts = counts.loc[counts.sum(axis=1) > 10]

tumor_cols  = metadata.loc[metadata["condition"] == "tumor",  "sample_id"].tolist()
normal_cols = metadata.loc[metadata["condition"] == "normal", "sample_id"].tolist()

# run t-tests across all genes — no seed set
pvals = {}
for gene in counts.index:
    t, p = stats.ttest_ind(
        counts.loc[gene, tumor_cols],
        counts.loc[gene, normal_cols]
    )
    pvals[gene] = p

pval_series = pd.Series(pvals)
sig_genes   = counts.loc[pval_series < 0.05]

# manually picked after inspecting the results
candidates = ["GENE_42", "GENE_107", "GENE_889", "GENE_1204"]

sig_genes.to_csv("significant_results.csv")
print(f"Done. {len(sig_genes)} significant genes found.")
```

> **Also check the repo ([RNA_seq_analysis](https://github.com/mitreacristina/RNAseq_analysis)):** No `requirements.txt`, no `environment.yml`, no `Dockerfile`, no `README`, no `.gitignore`. The only commit message is "final version". Raw FASTQ files are tracked directly in git. There are no unit tests.

---

#### Option B — R

```r
# analysis.R — DE analysis script
# (no author, no date, no version info)

setwd("/Users/mitreacristina/Desktop/rnaseq_project")

counts   <- read.csv("counts_final_FINAL_v3.csv")
metadata <- read.csv("meta.csv")

# filter low counts (threshold chosen after looking at results)
counts <- counts[rowSums(counts) > 10, ]

# run t-tests on all genes
pvals <- apply(counts, 1, function(x) {
  t.test(x[metadata$condition == "tumor"],
         x[metadata$condition == "normal"])$p.value
})

sig_genes <- counts[pvals < 0.05, ]

# pick top candidates by hand
candidates <- c("GENE_42", "GENE_107", "GENE_889", "GENE_1204")

write.csv(sig_genes, "significant_results.csv")
```

> **Also check the repo ([RNA_seq_analysis](https://github.com/mitreacristina/RNAseq_analysis)):** No `renv.lock` (note: `renv` supersedes the older `packrat` for R dependency management), no `Dockerfile` or Apptainer (formerly Singularity) definition, no `README`, no `.gitignore`. The only commit message is "final version". Raw FASTQ files are tracked directly in git. There are no unit tests. The required packages must be installed manually with no version guidance.

---

#### Task 2A — Code audit table

Complete the table for your chosen language.

| Issue Category | What's Wrong Here | How to Fix It |
|---|---|---|
| Environment / dependency management | | |
| Containerization | | |
| Hardcoded file paths | | |
| File naming / version control | | |
| Large files in git | | |
| Statistical test choice for count data | | |
| Multiple testing correction | | |
| Post-hoc filtering threshold | | |
| Manual candidate selection | | |
| Random seed / stochastic reproducibility | | |
| Code documentation | | |
| Unit testing | | |
| Commit message quality | | |

#### Task 2B — Rewrite the stats section

Rewrite the testing section (the loop or `apply` block) using a statistically appropriate method for RNA-seq count data.

RNA-seq counts are discrete, non-negative integers that exhibit **overdispersion** — the variance exceeds the mean, which violates both the normality assumption of the t-test and the equidispersion assumption of a Poisson model. The **negative binomial distribution** explicitly models this overdispersion and is the standard choice for RNA-seq differential expression.

Your rewrite must:

1. Use the correct package and key function calls with their main arguments — pseudocode is acceptable but must show the analysis pipeline steps, not just name the package
   - Python: `pydeseq2` (`DeseqDataSet` → `deseq2()` → `DeseqStats` → `results_df`)
   - R: `DESeq2` (`DESeqDataSetFromMatrix` → `DESeq` → `results`)
2. Apply Benjamini–Hochberg FDR correction and name the function used
3. Set a random seed (Python) or document why one is not needed (R)
4. Include at least one comment per major step explaining *why*, not just *what*

#### Task 2C — Environment reproducibility

Rather than simply writing an environment file, answer the following:

1. Name two specific things that would **silently break** if a collaborator ran the original script two years from now without any pinned versions. Be concrete — name a package and a type of change that could occur.
2. Write a minimal `environment.yml` (Python/Conda) or `renv` initialization (R) pinning at least 5 relevant packages to specific versions.
3. In one sentence: when is a pinned environment file alone *not* sufficient for full reproducibility, and what additional tool addresses this?

---

### Step 3 — FAIR & Bias Assessment *(~9 min)*

#### Task 3A — FAIR score

Rate each dimension on a 1–5 scale (1 = completely absent, 5 = fully compliant) and justify in one sentence. Refer to specific evidence from the methods excerpt.

| Dimension | Score (1–5) | Justification (cite specific evidence from the methods) |
|---|---|---|
| Findable | | |
| Accessible | | |
| Interoperable | | |
| Reusable | | |

> Consider: Does "available upon reasonable request" satisfy Accessibility? Does the absence of a `LICENSE` file affect Reusability? Would depositing raw data in NCBI Gene Expression Omnibus (GEO) or Zenodo change the Findability score?

#### Task 3B — Bias identification

List at least **three distinct sources of bias or confounding** present in this study. For each one provide:
1. The **type of bias** (e.g., selection bias, confirmation bias, technical confounding)
2. The **specific claim in the methods excerpt** that is undermined — do not write generic answers like "it threatens validity"
3. One concrete **mitigation strategy** with a named tool or approach

#### Task 3C — Team composition

The methods note that both authors are from the same lab and share the same clinical specialty (oncology). Answer both parts:

1. **Confirmation bias:** How does same-lab authorship increase the risk of confirmation bias, and what is one structural safeguard that can be put in place at the analysis stage?
2. **Domain blind spots:** What type of expertise is missing from this two-person oncology team that would have caught the core statistical flaw in this script, and how would you recruit that expertise at the study design phase?

---

### Step 4 — Power Calculation & Design Proposal *(~8 min)*

------------------------------
#### Task 4A — Was the study powered?

Cohort A has n = 6 per group. Use the formula below to estimate whether this study was adequately powered to detect a 2-fold change in expression.

> **Important caveat:** The formula below applies to a two-sample t-test on continuous, normally distributed data — the same (incorrect) test used in the script. This gives a lower-bound estimate. For a properly designed RNA-seq study, use purpose-built tools such as `RnaSeqSampleSize` (R/Bioconductor) or `RNASeqPower` (R), which use negative binomial models. The rule of thumb for pilot RNA-seq studies is a minimum of 3 biological replicates per group; for discovery studies, ≥10–20 per group is typical after accounting for multiple testing.

---

##### Part 1 — Worked example: t-test power (microarray context)

The t-test power formula is appropriate for **continuous, normally distributed data** — for example, log₂ microarray intensities. We work through it here to build intuition before moving to the correct model for RNA-seq.

You are designing a microarray study comparing tumor vs. normal tissue. From pilot data:

| Parameter | Symbol | Value |
|---|---|---|
| Significance threshold | α | 0.05 (two-tailed) |
| Desired power | 1 − β | 0.80 |
| Standard deviation of log₂ intensities | σ | 0.8 |
| Minimum difference to detect | δ | 1.0 (a 2-fold change on the log₂ scale) |

**Formula:**

$$n \approx 2 \times \left[\frac{(z_{\alpha/2} + z_\beta) \times \sigma}{\delta}\right]^2$$

where z_α/2 = 1.96 and z_β = 0.84.

**Calculation:**

$$n \approx 2 \times \left[\frac{(1.96 + 0.84) \times 0.8}{1.0}\right]^2 = 2 \times [2.24]^2 \approx \mathbf{11 \text{ per group}}$$

Verify with code — both snippets should return n = 11:

**R**
```r
library(pwr)
result <- pwr.t.test(
  d           = 1.0 / 0.8,   # Cohen's d = delta / sigma
  sig.level   = 0.05,
  power       = 0.80,
  type        = "two.sample",
  alternative = "two.sided"
)
ceiling(result$n)
```

**Python**
```python
import math
from statsmodels.stats.power import TTestIndPower

n = TTestIndPower().solve_power(
    effect_size = 1.0 / 0.8,  # Cohen's d = delta / sigma
    alpha       = 0.05,
    power       = 0.80,
    alternative = "two-sided"
)
print(math.ceil(n))
```

> **Why this formula does not apply to RNA-seq:** RNA-seq counts are discrete, non-negative integers whose variance scales with the mean — they follow a **negative binomial distribution** with Var(X) = μ + μ²·φ, where φ is the gene-specific dispersion (φ = BCV², the squared biological coefficient of variation). The t-test formula assumes constant variance and normality. More importantly, this is a high-throughput experiment testing 18,000 genes simultaneously — the single-gene formula has no concept of multiple testing and will severely underestimate the required sample size.

---

#### Part 2 — Correct power calculation for RNA-seq

Because we are selecting DE genes across the full transcriptome, the power calculation must account for the FDR correction applied across all 18,000 tests. Use `RnaSeqSampleSize` (R/Bioconductor), which models the negative binomial distribution and derives the per-gene significance threshold directly from the genome-wide FDR target (Zhao et al., *BMC Bioinformatics*, 2018).

**Parameters:**

| Parameter | Value | Notes |
|---|---|---|
| Total genes tested | 18,000 | As stated in the methods |
| Expected DE genes | 900 | 5% of 18,000 — a reasonable prior for tumor vs. normal |
| Fold change to detect | 2 | Minimum biologically meaningful change |
| Mean count (control) | 5 | Average read depth for a DE gene of interest |
| Dispersion (φ = BCV²) | 0.16 | BCV = 0.4, typical for human RNA-seq; φ = 0.4² = 0.16 |
| FDR threshold | 0.05 | Benjamini–Hochberg |
| Target power | 0.80 | |

**R**
```r
# install if needed: BiocManager::install("RnaSeqSampleSize")
library(RnaSeqSampleSize)

n <- sample_size(
  power   = 0.80,   # desired power (1 - Type II error rate)
  m       = 18000,  # total number of genes tested
  m1      = 900,    # expected number of truly DE genes
  f       = 0.05,   # FDR threshold (Benjamini-Hochberg)
  rho     = 2,      # fold change to detect (Treatment / Control)
  lambda0 = 5,      # mean read count for a DE gene in the control group
  phi0    = 0.16,   # dispersion = BCV^2 = 0.4^2; typical for human RNA-seq
  k       = 1,      # ratio of group sizes (1 = balanced design)
  w       = 1       # ratio of normalization factors (1 = no systematic bias)
)
print(n)
```

**Python** — call the R function directly via `rpy2`:
```python
import math
from rpy2.robjects.packages import importr
from rpy2 import robjects

sample_size = importr("RnaSeqSampleSize").sample_size

n = sample_size(
    power   = robjects.FloatVector([0.80]),
    m       = robjects.IntVector([18000]),
    m1      = robjects.IntVector([900]),
    f       = robjects.FloatVector([0.05]),
    rho     = robjects.FloatVector([2]),
    lambda0 = robjects.FloatVector([5]),
    phi0    = robjects.FloatVector([0.16]),
    k       = robjects.FloatVector([1]),
    w       = robjects.FloatVector([1])
)
print(math.ceil(float(n[0])))
```

> **Expected result:** Approximately **20–25 samples per group** — consistent with published benchmarks for human RNA-seq studies with BCV = 0.4 and a 2-fold detection threshold. This is three to four times the n = 6 used in Cohort A, and five times the n = 4 in Cohort B.

**Answer these questions:**

1. What minimum n per group does `sample_size()` return with these parameters?
2. How does the answer change if you assume only 1% of genes are truly DE (m1 = 180)? What does this tell you about how sensitive the power estimate is to your assumptions?
3. The scenario uses n = 6 (Cohort A) and n = 4 (Cohort B). Beyond missing true DE genes, what is one additional scientific consequence of running an underpowered high-throughput study?


#### Task 4B — Redesign in 5 bullets

Write exactly **5 bullet points** describing an improved study design, one per course principle listed below. Each bullet must include: (a) one named tool or approach, (b) one sentence explaining which specific flaw it addresses and why, and (c) one limitation — something this fix does *not* solve.

Principles to address (one each): **statistics · software engineering · FAIR · bias/confounding · team diversity**

> Python tools to draw from: Snakemake, Zenodo, Docker, GEO, `pydeseq2`, `pytest`, DVC, OSF  
> R tools to draw from: Nextflow, Zenodo, Apptainer, GEO, `DESeq2`, `testthat`, `renv`, `RnaSeqSampleSize`

---

### Step 5 — Synthesis *(~7 min)*

#### Task 5A — Reflection

In 3–5 sentences: *Which single flaw in this study do you think is most common in published bioinformatics literature, and why is it so persistent?* Draw on at least one concept from the course. There is no single correct answer — we are looking for critical reasoning and specificity.

#### Task 5B — Prioritization *(~7 min)*

You have one week before your collaborator submits the grant. You cannot fix everything. **Choose the single most important flaw to address first.** In 3–5 sentences, defend your choice. Your answer must: (1) name the flaw, (2) explain what harm it causes if left unfixed, and (3) acknowledge the strongest counter-argument for fixing a different flaw instead.

> **Optional extension:** Find one real published paper and identify one reproducibility or FAIR issue in its methods or data availability statement. Write two sentences. 

---

### Deliverables *(~10 min for exercise prep and deliverables check)*

Submit all of the following as a pull request to this repository or single PDF or `.md` file via the course portal.

1. Completed annotation checklist (Task 1A) with one sentence per item
2. Completed code audit table (Task 2A), rewritten stats section (Task 2B), and environment reproducibility answers (Task 2C)
3. FAIR scores with evidence-based justifications (Task 3A), bias analysis (Task 3B), and team composition answers (Task 3C parts 1 and 2)
4. Power calculation with workthrough shown and genome-wide interpretation (Task 4A), plus 5-bullet redesign with tool, rationale, and limitation per bullet (Task 4B)
5. Reflection (Task 5A) and prioritization argument (Task 5B)
