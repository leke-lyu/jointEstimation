# Discrete Trait Analysis — County

This module performs Bayesian phylogeographic discrete trait analysis (DTA) to infer transmission rates between the 9 counties of Greater Houston, then computes epidemiological metrics (Source Sink Score, Local Import Score, Persistence Time) from the posterior trees.

For full methodological details, see [Lyu et al. (2025)](https://doi.org/10.1073/pnas.2501435122).

## Pipeline overview

```
00_prepareData.sh   Symlink cluster data from introductionAnalysis
        ↓
01_generateXML.sh   Generate BEAST XML with county-level traits
        ↓
02_runBeast.sh      Submit 5 independent BEAST runs (SLURM)
        ↓
03_logCombine.sh    Combine log files AND tree files with logcombiner
        ↓
04_divideAndConquer.sh   Split combined tree files → compute metrics per posterior tree
        ↓
05_localPhylogeographicAnalysis.sh   Aggregate metrics → generate map and metric plots
```

## Geographic scope

The 9 counties of Greater Houston:

| County | Role in paper |
|--------|--------------|
| Harris | Primary viral source (SSS ~ 0.93) |
| Fort Bend | Top recipient from Harris |
| Montgomery | Second-largest recipient from Harris |
| Brazoria | Recipient from Harris |
| Galveston | Intermediate |
| Liberty | Peripheral |
| Waller | Peripheral |
| Chambers | Peripheral |
| Austin | Peripheral (smallest sample) |

## Prerequisites

Outputs from the `introductionAnalysis` stage:
- `clusters.RData` — extracted cluster trees
- `filteredHoustonDeltaMeta.RData` — metadata with county assignments

---

## 00_prepareData.sh

Creates symlinks to shared data files from `introductionAnalysis`.

---

## 01_generateXML.sh

### Objective

Generate the BEAST XML for county-level phylogeographic inference.

### What it does

Runs `scripts/xml.R`, which:

1. Loads cluster trees and metadata.
2. Maps each tip to its county trait (spaces replaced with underscores, e.g., `Fort_Bend_County`).
3. Writes per-cluster NEXUS tree files to `tre/`.
4. Generates `output2.xml` — a BEAST XML that jointly estimates transition rates across all 82 clusters with:
   - **Non-reversible CTMC** with 72 rate parameters (9 counties x 8 transitions each)
   - **BSSVS** with binary indicator variables
   - **Poisson prior** (mean=1, offset=0)

### BEAST run parameters

| Parameter | Value |
|-----------|-------|
| Chain length | 500,000,000 (longer due to more parameters) |
| Log every | 50,000 states |
| Empirical trees per cluster | 500 |

### Outputs

| File | Description |
|------|-------------|
| `output2.xml` | BEAST XML configuration file |
| `tre/*.nex` | Per-cluster NEXUS tree files |

---

## 02_runBeast.sh

Submits 5 independent BEAST runs via SLURM (same setup as DTAAge: 2 CPUs, 2 GB, 720-hour limit).

---

## 03_logCombine.sh

### Objective

Combine posterior samples from all 5 runs, for both log files and tree files.

### What it does

Uses `logcombiner` with:
- **Burn-in:** 100,000,000 states
- **Resample:** every 1,600,000 states
- Combines **log files** into `output/combined/combined_*.log`
- Combines **tree files** (`.trees`) into `output/combined/combined_*.trees` — these contain the posterior trees with inferred ancestral county assignments, needed for computing transmission metrics.

---

## 04_divideAndConquer.sh

### Objective

Split the combined posterior tree files into individual trees and compute Source Sink Score, Local Import Score, and Persistence Time for each posterior sample.

### What it does

1. For each combined `.trees` file, uses `awk` to extract individual tree states into separate files: `{i}_{cladeName}.tree` (up to 1,000 posterior samples).
2. For each posterior sample index `i`, submits a SLURM job running `scripts/processTrees.R` on all cluster trees at that index.

### scripts/processTrees.R

For each posterior tree:

1. Reads the annotated tree and extracts parent–child `location` assignments.
2. **Persistence Time:** For each tip, traces backward through the tree until the ancestral location changes. The persistence time is the cumulative branch length from the tip to the point where the lineage entered the current location (using half the transition branch). Result is in years.
3. **Source Sink Score (SSS):**
   ```
   SSS = (exports - imports) / (exports + imports)
   ```
   - Imports = branches where `from != county` and `to == county`
   - Exports = branches where `from == county` and `to != county`
4. **Local Import Score (LIS):**
   ```
   LIS = imports / (imports + local_transmission)
   ```
   - Local transmission = branches where `from == county` and `to == county`

5. Saves a data frame with `Location`, `MedianLastingTime`, `sss`, `lis` per posterior sample.

### Outputs

| File | Description |
|------|-------------|
| `output/combined/{i}.RData` | Per-posterior-sample metrics for all counties (1,000 files) |

---

## 05_localPhylogeographicAnalysis.sh

### Objective

Aggregate the per-posterior-sample metrics and generate publication figures showing county-level transmission patterns.

### What it does

Runs `scripts/localPhylogeographicAnalysis.R`, which:

1. Loads all `.RData` files from `output/combined/` and aggregates across posterior samples.
2. For each county, computes the **median, min, and max** of SSS, LIS, and Persistence Time across the posterior distribution. These ranges represent the 95% HPD interval.
3. Reads the combined BEAST log to compute **Bayes Factors** for all pairwise county transitions (same method as DTAAge).
4. Filters transitions with BF > 100 (decisive support) for the geographic map.
5. Generates two figures.

### Interpreting the results

**`map.pdf` — Geographic diffusion map:**
- Shows the 9 Greater Houston counties as blue polygons.
- **Arrows** indicate supported transmission pathways (BF > 100 only).
- **Arrow thickness** is proportional to the estimated transition rate.
- **Key finding:** Harris County is the dominant source, with strong outward arrows to Fort Bend, Montgomery, and Brazoria. Few arrows point inward to Harris, confirming its role as the primary viral hub.

**`metric.pdf` — Epidemiological metrics (3 stacked panels):**

Counties are ordered on the X-axis: Harris, Fort Bend, Montgomery, Brazoria, Galveston, Liberty, Waller, Chambers, Austin.

**Top panel — Source Sink Score:**
- **Y-axis:** SSS (−1 to +1)
- **Green diamonds:** Median SSS per county
- **Golden error bars:** Range across posterior samples (95% HPD)
- Harris has a strongly positive SSS (~0.93), confirming it as the dominant viral source. Surrounding counties have negative scores, indicating they are net sinks.

**Middle panel — Local Import Score:**
- **Y-axis:** LIS (0 to 1)
- **Orange squares:** Median LIS per county
- **Golden error bars:** 95% HPD range
- Harris has a near-zero LIS (~0.004), meaning its epidemic was almost entirely sustained by local transmission. Peripheral counties have higher LIS values, indicating greater reliance on continued introductions from Harris.

**Bottom panel — Persistence Time:**
- **Y-axis:** Median persistence time (years, 0 to ~0.6)
- **Blue circles:** Median persistence time per county
- **Golden error bars:** 95% HPD range
- Harris has the longest persistence time (~0.59 years), meaning viral lineages circulated locally for extended periods. Peripheral counties show shorter persistence, consistent with repeated introductions that do not establish sustained local chains.
- **Key insight from the paper:** Regions with high SSS (sources) tend to have low LIS and long persistence times. A well-established, self-sustaining local epidemic is what enables a region to function as a viral source.

**`output_tb.csv`:**
- Bayes Factor data frame with geographic coordinates, suitable for custom mapping.

### Outputs

| File | Description |
|------|-------------|
| `metric.pdf` | 3-panel bar chart (SSS, LIS, Persistence Time) |
| `map.pdf` | Geographic map with directed transmission arrows |
| `output_tb.csv` | BF data frame with coordinates |

---

## R dependencies

| Package | Used by |
|---------|---------|
| `treeio`, `ggtree` | processTrees.R |
| `dplyr`, `magrittr`, `stringr` | All scripts |
| `sf`, `tigris` | localPhylogeographicAnalysis.R (county boundaries) |
| `ggplot2`, `ggpubr` | localPhylogeographicAnalysis.R |
| `lubridate` | xml.R |
| `ape` | xml.R |
