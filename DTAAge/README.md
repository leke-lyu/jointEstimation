# Discrete Trait Analysis — Age Groups

This module performs Bayesian phylogeographic discrete trait analysis (DTA) to infer transmission rates between age groups across all 82 introduction clusters, using a joint estimation framework in BEAST v1.10.4.

For full methodological details, see [Lyu et al. (2025)](https://doi.org/10.1073/pnas.2501435122).

## Pipeline overview

```
00_prepareData.sh   Symlink cluster data from introductionAnalysis
        ↓
01_generateXML.sh   Generate BEAST XML with age group traits
        ↓
02_runBeast.sh      Submit 5 independent BEAST runs (SLURM)
        ↓
03_logCombine.sh    Combine MCMC outputs with logcombiner
        ↓
04_plot.sh          Compute Bayes Factors and generate chord diagram
```

## Prerequisites

This module requires outputs from the `introductionAnalysis` stage:
- `clusters.RData` — extracted cluster trees
- `filteredHoustonDeltaMeta.RData` — metadata with age groups

---

## 00_prepareData.sh

### Objective

Create symlinks to the shared data files from `introductionAnalysis`.

### What it does

Creates a `data/` directory and symlinks `clusters.RData` and `filteredHoustonDeltaMeta.RData` from the introduction analysis stage.

---

## 01_generateXML.sh

### Objective

Generate the BEAST XML configuration file and per-cluster NEXUS tree files for the age group DTA.

### What it does

Runs `scripts/xml.R`, which:

1. Loads cluster trees and metadata.
2. Maps each tip to its age group trait (spaces replaced with underscores: `Infants_and_Children`, `Teenagers`, `Young_Adults`, `Middle-aged_Adults`, `Seniors`).
3. Writes each cluster's phylogenetic tree as a separate NEXUS file in `tre/`.
4. Generates a single BEAST XML (`output.xml`) that jointly estimates transition rates across all 82 clusters. The XML includes:
   - **Taxa blocks** with decimal sampling dates
   - **Alignment blocks** with dummy sequences (NNNN) — the analysis is trait-only, not sequence-based
   - **Empirical tree distributions** — one per cluster, using 500 sampled trees
   - **Discrete trait model**: Non-reversible continuous-time Markov chain (CTMC) with 20 rate parameters (5 age groups x 4 possible transitions each)
   - **Bayesian stochastic search variable selection (BSSVS)** with binary indicator variables to identify which transitions have non-zero rates
   - **Prior**: Poisson prior (mean=1, offset=0) on the total number of non-zero rates
   - **MCMC operators**: scale, bit-flip, delta-exchange

### BEAST run parameters

| Parameter | Value |
|-----------|-------|
| Chain length | 100,000,000 |
| Log every | 10,000 states |
| Empirical trees per cluster | 500 |

### Outputs

| File | Description |
|------|-------------|
| `output.xml` | BEAST XML configuration file |
| `tre/*.nex` | Per-cluster NEXUS tree files |

---

## 02_runBeast.sh

### Objective

Submit 5 independent BEAST MCMC runs for convergence assessment.

### What it does

For each XML file found in the directory:
1. Creates a subdirectory with 5 run folders (`run1`–`run5`).
2. Generates and submits a SLURM job script for each run.

### Computational requirements

| Parameter | Value |
|-----------|-------|
| Partition | bahl_p |
| CPUs per run | 2 |
| Memory | 2 GB |
| Time limit | 720 hours (30 days) |
| Software | BEAST v1.10.4 |

---

## 03_logCombine.sh

### Objective

Combine the posterior samples from all 5 independent BEAST runs into a single log file after discarding burn-in.

### What it does

Uses BEAST's `logcombiner` tool to:
1. Discard the first 20,000,000 states as burn-in from each run.
2. Resample every 40,000 states.
3. Combine log files from all 5 runs into `output/combined/combined_All_clades.log`.

---

## 04_plot.sh

### Objective

Compute Bayes Factors for all pairwise age group transitions and visualize supported transmission pathways as a chord diagram.

### What it does

Runs `scripts/BF.R`, which:

1. Reads `output/combined/combined_All_clades.log`.
2. Extracts BSSVS indicator and rate columns for each pairwise transition.
3. Computes **actual transition rates** = raw rate x indicator (zero when indicator is off).
4. Computes **Bayes Factors** for each transition using:
   ```
   BF = [p * (1 - q)] / [(1 - p) * q]
   ```
   where `p` = posterior mean of the indicator, and `q` = prior expectation = `(log(2) + k - 1) / (k * (k - 1))` with `k` = number of age groups (5).
5. Filters transitions with posterior mean of indicator > 0.5.
6. Categorizes Bayes Factor support:

| Category | BF threshold | Interpretation |
|----------|-------------|----------------|
| type_a | > 100 | Decisive support |
| type_b | > 30 | Very strong support |
| type_c | > 10 | Strong support |
| type_d | > 3 | Positive support |

7. Generates a **chord diagram** (`flow.pdf`) showing directional transmission flows between age groups.

### Interpreting the results

**`flow.pdf` — Age group transmission chord diagram:**
- **Sectors** around the circle represent age groups, split into source (left, `src-`) and sink (right, `snk-`) sides.
- **Chords** connect source age groups to sink age groups. Chord thickness is proportional to the estimated transition rate.
- **Chord color** indicates the strength of Bayes Factor support (darker = stronger evidence).
- **Key finding from the paper:** Young adults and middle-aged adults are the primary drivers of transmission. The strongest supported transitions (BF > 100) include young → middle-aged, middle-aged → young, middle-aged → children, middle-aged → seniors, and middle-aged → teenagers. Seniors primarily receive virus from younger age groups.

**`bf.RData`:**
- Contains the `bf_df` data frame with columns: `from`, `to`, `posterior.mean.of.indicators`, `posterior.mean.of.location.rates`, `TransitionRate`, `bf`, `group`.
- Load with `load("bf.RData")` for custom analysis.

### Outputs

| File | Description |
|------|-------------|
| `bf.RData` | Bayes Factor data frame for all transitions |
| `flow.pdf` | Chord diagram of age group transmission pathways |

---

## R dependencies

| Package | Used by |
|---------|---------|
| `dplyr`, `magrittr`, `stringr` | xml.R, BF.R |
| `lubridate` | xml.R |
| `ape` | xml.R |
| `ggplot2` | BF.R |
| `circlize` | BF.R (chord diagram) |
