# Tip-Trait Association Testing (TTAT)

This module tests whether demographic traits (age group, sex) are non-randomly distributed across the phylogenetic trees of the 82 introduction clusters. A significant association indicates that viral transmission is more constrained within certain trait groups than expected by chance — i.e., that population structure influences disease spread.

For full methodological details, see [Lyu et al. (2025)](https://doi.org/10.1073/pnas.2501435122).

The TTAT method is implemented as an R package: [github.com/leke-lyu/TTAT](https://github.com/leke-lyu/TTAT).

## Pipeline overview

```
00_prepareData.sh   Symlink cluster data from introductionAnalysis
        ↓
01_ttat.sh          Run trait-tree association tests for age groups and sex
```

## Prerequisites

Outputs from the `introductionAnalysis` stage:
- `clusters.RData` — extracted cluster trees
- `filteredHoustonDeltaMeta.RData` — metadata with age groups and sex

---

## 00_prepareData.sh

Creates symlinks to shared data files from `introductionAnalysis`.

---

## 01_ttat.sh

### Objective

Quantify the degree to which age group and sex are phylogenetically clustered within each of the 82 introduction clusters, and assess statistical significance via permutation testing.

### What it does

Runs `scripts/ttat.R`, which:

1. Loads cluster trees and metadata.
2. For each trait (age groups, then sex), iterates over all 82 clusters:
   - Drops tips with `"Unknown"` trait values.
   - Skips clusters with fewer than 2 unique trait values or fewer than 10 remaining tips.
   - Computes the **Association Index (AI)** using the `TTAT` package:
     - AI measures the degree of phylogenetic clustering of trait values. Lower AI = stronger clustering (traits are more segregated on the tree).
   - Generates a **null distribution** by permuting trait labels 1,000 times while preserving the tree topology and trait frequency distribution.
   - Computes the **p-value**: proportion of null AI values <= observed AI.
3. Saves results and generates visualizations.

### How the Association Index works

The Association Index quantifies whether tips with the same trait tend to cluster together on the phylogeny more than expected by chance. For each internal node, it measures how evenly trait categories are distributed among its descendants. The overall AI is a weighted sum across all internal nodes.

- **Low AI** → tips with similar traits cluster together (non-random; population structure influences transmission)
- **High AI** → traits are randomly distributed across the tree (no evidence of structured transmission)

### Interpreting the results

**`ttat.pdf` — Scatter plot of all 82 clusters:**
- **X-axis:** Cluster size (log scale)
- **Y-axis:** P-value from the permutation test
- **Blue circles:** Age group tests
- **Green triangles:** Sex tests
- **Red dashed line:** p = 0.05 significance threshold
- Points **below** the red line indicate clusters where the trait is significantly clustered on the phylogeny.
- **Key finding:** 20 of 82 clusters showed significant age group associations (p < 0.05), but only 6 showed significant sex associations. This indicates that "viral transmission is more constrained within age groups than sex groups," consistent with age-structured social environments (schools, workplaces, nursing homes) shaping transmission more than sex-based mixing patterns.

**`nulldistribution.pdf` — Example null distribution (cluster 24):**
- **X-axis:** Association Index values
- **Y-axis:** Frequency
- **Blue histogram:** Distribution of AI values from 1,000 random permutations
- **Red dashed line:** Observed AI value for this cluster
- When the observed AI falls far to the left of the null distribution, the p-value is small — the trait is more clustered than expected by chance.

**`mycluster.pdf` — Example phylogenetic tree (cluster 24):**
- Time-scaled phylogenetic tree with tip colors indicating age group.
- Visually demonstrates age group clustering: tips of similar colors tend to group together on branches, reflecting within-age-group transmission chains.
- Color scheme: blue = Young Adults, pink = Middle-aged Adults, orange = Infants and Children, purple = Seniors, light pink = Teenagers.

**`ttat.RData`:**
- Contains two data frames:
  - `age_groups`: cluster_size and p_value for each cluster (age group test)
  - `sex`: cluster_size and p_value for each cluster (sex test)
- Load with `load("ttat.RData")` for custom analysis.

### Outputs

| File | Description |
|------|-------------|
| `ttat.RData` | P-values for age group and sex tests across all 82 clusters |
| `ttat.pdf` | Scatter plot of cluster size vs. p-value |
| `nulldistribution.pdf` | Null distribution histogram for example cluster 24 |
| `mycluster.pdf` | Phylogenetic tree for cluster 24 colored by age group |

---

## R dependencies

| Package | Used by |
|---------|---------|
| `TTAT` | ttat.R (install from [github.com/leke-lyu/TTAT](https://github.com/leke-lyu/TTAT)) |
| `dplyr`, `magrittr` | ttat.R |
| `ggplot2` | ttat.R |
| `ggtree` | ttat.R (cluster tree visualization) |
| `ape` | ttat.R (tree manipulation) |
