# PhyloDeep — Phylodynamic Parameter Estimation

This module uses [PhyloDeep](https://github.com/evolbioinfo/phylodeep), a deep learning-based phylodynamic inference tool, to estimate the basic reproduction number (R0) and infectious period for each introduction cluster from its tree topology.

## Pipeline overview

```
run.sh               Extract Newick trees → run paramdeep on each cluster
        ↓
summarize_clades.R   Compile R0 and infectious period into a summary table
```

## Prerequisites

- Cluster tree files (`.tre` format) in the `tre/` directory. These are the NEXUS trees from `introductionAnalysis/clusters.RData`, exported as individual Newick tree files.
- **Conda environment** `phyloenv` with `paramdeep` installed.

---

## run.sh

### Objective

Estimate phylodynamic parameters for each introduction cluster using the Birth-Death model with a feedforward neural network (BD-FFNN).

### What it does

1. Activates the `phyloenv` conda environment.
2. For each `.tre` file in `tre/`:
   - Extracts the Newick string from the NEXUS-formatted tree file (strips the NEXUS header, keeping only the tree topology after `[&R]`).
   - Overwrites the file with the clean Newick string.
   - Runs `paramdeep` with parameters:
     - `-t`: Input tree file
     - `-p 0.1`: Sampling probability (fraction of infected individuals that were sampled)
     - `-m BD`: Birth-Death epidemiological model
     - `-v FFNN_SUMSTATS`: Use feedforward neural network trained on summary statistics
   - Saves output to `{clade}_BD_FFNN.csv`.

### paramdeep model

PhyloDeep uses pre-trained neural networks to infer epidemiological parameters directly from phylogenetic tree summary statistics, bypassing the need for computationally expensive Bayesian MCMC inference. The Birth-Death (BD) model parameterizes the epidemic process with:
- **R0 (basic reproduction number):** Average number of secondary infections per primary case
- **Infectious period:** Average duration of infectiousness

### Outputs

| File | Description |
|------|-------------|
| `tre/{Clade}_BD_FFNN.csv` | Per-cluster estimates of R0 and infectious period |

---

## summarize_clades.R

### Objective

Compile the per-cluster PhyloDeep results into a single summary table.

### What it does

1. Finds all `*_BD_FFNN.csv` files in `tre/`.
2. For each cluster:
   - Reads the corresponding `.tre` file to count the number of tips (cluster size).
   - Extracts R0 and infectious period from the CSV.
3. Produces a summary table with one row per cluster.

### Interpreting the results

**`clade_summary.csv`:**

| Column | Description |
|--------|-------------|
| `Clade` | Cluster identifier (e.g., `Clade_17`) |
| `Size` | Number of tips (isolates) in the cluster |
| `R_naught` | Estimated basic reproduction number. Values > 1 indicate a growing epidemic within the cluster; values < 1 indicate a declining one |
| `Infectious_period` | Estimated average duration of infectiousness (in the time units of the tree, typically years) |

- Larger clusters that established sustained local transmission are expected to have R0 > 1.
- The infectious period estimates can be compared across clusters to assess whether transmission dynamics varied over the course of the Delta wave.

### Outputs

| File | Description |
|------|-------------|
| `tre/clade_summary.csv` | Summary table: Clade, Size, R_naught, Infectious_period |

---

## Software dependencies

| Tool | Source |
|------|--------|
| `paramdeep` (PhyloDeep) | [github.com/evolbioinfo/phylodeep](https://github.com/evolbioinfo/phylodeep) |
| Conda environment | `phyloenv` |
| R packages: `ape`, `dplyr` | summarize_clades.R |
