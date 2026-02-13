# Introduction Analysis

This is the foundational stage of the pipeline. It processes a posterior set of phylogeographic trees from BEAST to identify viral introduction events into Greater Houston, extract locally circulating clusters, and prepare metadata for all downstream analyses.

For full methodological details, see [Lyu et al. (2025)](https://doi.org/10.1073/pnas.2501435122).

## Pipeline overview

```
01_divideAndConquer.sh   Split combined.trees into individual trees → summarize introductions per tree
        ↓
02_importationAnalysis.sh   Aggregate introduction summaries → plots and statistics
        ↓
03_extractClusters.sh   Extract clusters from representative tree (tree137)
        ↓
04_prepareMetaData.sh   Map metadata (age, sex, county) to cluster tips
        ↓
05_plotClusters.sh   Visualize selected clusters colored by county
```

All downstream modules (DTAAge, DTACounty, GLM, TTAT, phylodeep) depend on the outputs of this stage.

---

## 01_divideAndConquer.sh

### Objective

Split the combined posterior tree file into individual NEXUS trees and summarize introduction events for each, enabling posterior uncertainty quantification.

### What it does

1. Reads `data/combined.trees` — a file containing multiple posterior trees from a BEAST phylogeographic run.
2. Uses `awk` to extract each tree state into a separate NEXUS file (`tree/tree1.nexus`, `tree/tree2.nexus`, ...).
3. For each extracted tree, submits a SLURM job that runs `scripts/importTreeSummary.R` with:
   - The tree file
   - Focal area: `"Texas"`
   - Maximum date: `2021.81369863014` (decimal date for the latest sample)

### scripts/importTreeSummary.R

For a given tree:
1. Reads the NEXUS tree and converts it to a data frame.
2. Identifies **introduction events** — branches where the parent node's `location` is not the focal area but the child node's `location` is. The introduction time is estimated as the midpoint of that branch.
3. For each introduction, extracts the downstream monophyletic subtree to determine the cluster size (number of Houston tips descended from that introduction).
4. Records the source of each introduction (`"Domestic"` or `"International"`), the epidemiological week, and the cluster size.
5. Saves per-tree results as `tree/{treeName}.RData`.

### Outputs

| File | Description |
|------|-------------|
| `tree/tree{N}.nexus` | Individual posterior trees extracted from combined.trees |
| `tree/tree{N}.RData` | Introduction summary per tree: node, source, epi-week, cluster size |

---

## 02_importationAnalysis.sh

### Objective

Aggregate the per-tree introduction summaries into a single data frame and generate publication figures characterizing the dynamics of viral introduction into Houston.

### What it does

Runs `scripts/posteriorSummary.R`, which:

1. Loads all `.RData` files from the `tree/` directory and combines them into a single data frame (`wholeDF`).
2. Categorizes each introduction by cluster size: singletons, clusters with size <= 10, clusters with size > 10.
3. Computes summary statistics across posterior trees:
   - Total number of introductions (median, range across trees)
   - Number of singletons
   - Number of domestic vs. international introductions
4. Generates two figures:
   - **`sizeSummary.pdf`**: Weekly introduction frequency broken down by cluster size category. Lines show medians; shaded ribbons show the range across posterior trees (95% HPD).
   - **`DvsI.pdf`**: Weekly introduction frequency colored by source (domestic vs. international). Same median/ribbon format.
5. Fits **negative binomial distributions** to cluster sizes (overall, domestic-only, international-only) to estimate the mean cluster size (r) and dispersion parameter (k) for each posterior tree.

### Interpreting the results

**`sizeSummary.pdf` — Introduction dynamics by cluster size:**
- **X-axis:** Date (2021)
- **Y-axis:** Number of introductions per epidemiological week
- **Grey line/ribbon:** Total imports (median and range)
- **Blue:** Singletons — introductions that did not establish onward transmission
- **Orange:** Small clusters (size <= 10)
- **Green:** Large clusters (size > 10)
- Early in the epidemic, a larger fraction of introductions seeded sizable clusters. Later (after the peak), most introductions resulted in singletons, reflecting growing population immunity.

**`DvsI.pdf` — Domestic vs. international introductions:**
- Domestic introductions dominate throughout, with international introductions peaking earlier. This indicates Houston was primarily seeded from other U.S. locations rather than serving as a direct international entry point for the Delta variant.

### Outputs

| File | Description |
|------|-------------|
| `wholeDF.RData` | Combined data frame of all introduction events across all posterior trees |
| `sizeSummary.pdf` | Weekly introductions by cluster size category |
| `DvsI.pdf` | Weekly introductions by domestic vs. international source |

---

## 03_extractClusters.sh

### Objective

From a representative posterior tree, identify introduction events and extract the locally circulating cluster subtrees for detailed downstream analysis.

### What it does

Runs `scripts/extractClusters.R`, which:

1. Reads `tree/tree137.nexus` — chosen as the representative tree (the tree with the median number of introductions).
2. Identifies introduction nodes: branches where `location` transitions from non-`"Texas"` (parent) to `"Texas"` (child).
3. For each introduction node, recursively extracts the downstream subtree containing all `"Texas"` descendants — this is the locally circulating cluster.
4. Filters clusters by size: keeps only those with **>10 tips** (isolates). This yields 82 clusters containing 6,455 sequences.
5. Cleans each cluster tree:
   - Removes nodes with no descendants
   - Collapses single-child internal nodes (adds their branch length to the child)
   - Re-indexes nodes for valid phylogenetic tree structure
6. Saves the results.

### Outputs

| File | Description |
|------|-------------|
| `clusters.RData` | Contains `newtreList` (list of 82 cleaned phylogenetic trees), `tipList` (tip labels per cluster), and `importSize` (cluster sizes) |

---

## 04_prepareMetaData.sh

### Objective

Link epidemiological metadata (age, sex, county) to the tips in the extracted clusters.

### What it does

Runs `scripts/prepareMetaData.R`, which:

1. Loads `clusters.RData` and `data/houstonDeltaMeta.RData` (raw metadata with GISAID names, dates, ZIP codes, age, sex).
2. Filters metadata to include only tips present in the 82 clusters.
3. Maps ZIP codes to Texas counties using the `zipcodeR` package.
4. Creates age group categories:

| Age group | Age range |
|-----------|-----------|
| Infants and Children | 0–12 |
| Teenagers | 13–19 |
| Young Adults | 20–35 |
| Middle-aged Adults | 36–55 |
| Seniors | 56+ |
| Unknown | Missing age data |

5. Prints frequency tables for age groups, sex, and county.

### Outputs

| File | Description |
|------|-------------|
| `filteredHoustonDeltaMeta.RData` | Metadata data frame with columns: GISAID_name, zip, collection_date, age, sex, county, new_age_groups, decimal_date |

---

## 05_plotClusters.sh

### Objective

Visualize selected cluster phylogenies colored by county to illustrate geographic mixing within clusters.

### What it does

Runs `scripts/plotClusters.R`, which:

1. Loads `clusters.RData` and `filteredHoustonDeltaMeta.RData`.
2. Defines a color scheme for the 9 Greater Houston counties:
   - Harris (blue), Fort Bend (pink), Montgomery (purple), Brazoria (orange), Galveston (coral), Liberty (teal), Waller (light pink), Chambers (dark blue), Austin (black)
3. Creates time-scaled phylogenetic trees for 5 selected clusters (78, 53, 45, 54, 10) with tips colored by county.
4. Stacks the 5 trees vertically in a combined figure.

### Interpreting the results

**`clusters.pdf`:**
- Each panel shows one cluster's phylogenetic tree.
- **X-axis:** Time (years, approximately March–October 2021)
- **Tip colors:** County of origin. A mix of colors within a cluster indicates cross-county transmission. Clusters dominated by a single color suggest geographically contained transmission.
- Larger clusters (top panels) tend to show more geographic mixing, while smaller clusters may be more geographically focused.

### Outputs

| File | Description |
|------|-------------|
| `clusters.pdf` | 5-panel phylogenetic tree visualization |
| `dummy_legend.pdf` | Standalone legend for county colors |

---

## scripts/bigtree.R

### Objective

A utility script for visualizing the full phylogeographic tree and marking introduction nodes. Useful for exploratory analysis.

### What it does

1. Reads `tree137.nexus` and identifies introduction nodes (same logic as `extractClusters.R`).
2. Plots the full tree with introduction nodes highlighted as yellow diamond markers.

---

## R dependencies

| Package | Used by |
|---------|---------|
| `treeio` | importTreeSummary.R, extractClusters.R, bigtree.R |
| `tidytree` | extractClusters.R |
| `ggtree` | plotClusters.R, bigtree.R |
| `dplyr`, `magrittr` | All scripts |
| `lubridate` | importTreeSummary.R |
| `ggplot2`, `ggpubr` | posteriorSummary.R |
| `MASS` | posteriorSummary.R (negative binomial fitting) |
| `zipcodeR` | prepareMetaData.R |
| `patchwork` | plotClusters.R |
| `ape` | extractClusters.R |

## Computational requirements

- **HPC:** SLURM scheduler (UGA cluster, `bahl_p` partition)
- **Modules:** `R/4.3.1-foss-2022a`
- **Per-tree jobs:** 1 CPU, 1 GB memory, 48-hour time limit
