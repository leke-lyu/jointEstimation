# GLM Analysis

This module extends the county-level discrete trait analysis with a Generalized Linear Model (GLM) to identify epidemiological and socioeconomic predictors that drive viral diffusion between counties. The GLM is embedded within BEAST's phylogeographic framework, where viral diffusion rates among counties serve as the outcome variable.

For full methodological details, see [Lyu et al. (2025)](https://doi.org/10.1073/pnas.2501435122).

## Pipeline overview

```
pca.R                  PCA-based variable selection to handle collinearity
        ↓
xml.R                  Generate BEAST XML with GLM substitution model
        ↓
02_runBeast.sh         Submit 5 independent BEAST runs (SLURM)
        ↓
03_logCombine.sh       Combine MCMC outputs with logcombiner
        ↓
04_glmSummary.sh       Summarize GLM coefficients, HPD intervals, and Bayes Factors
```

## Prerequisites

Outputs from the `introductionAnalysis` stage:
- `clusters.RData` — extracted cluster trees
- `filteredHoustonDeltaMeta.RData` — metadata with county assignments

---

## pca.R

### Objective

Perform variable selection via Principal Component Analysis (PCA) to address collinearity among candidate predictors before including them in the GLM.

### What it does

1. Defines 6 candidate predictor variables for the 9 Greater Houston counties:
   - **Cluster variables** (correlated with each other): case count, death count, population size, sample count
   - **Fixed variables** (always retained): median household income, population density
2. Performs PCA on all variables to visualize collinearity structure.
3. Scores each cluster variable by its maximum absolute correlation with the fixed variables.
4. Selects the cluster variable with the **lowest** correlation to the fixed variables — this minimizes collinearity in the final predictor set.

### Interpreting the results

The script outputs a ranking of cluster variables by collinearity with the fixed variables. The recommended variable (lowest max |r|) is selected to join median household income and population density in the final GLM predictor set. In practice, sample size was selected as the representative cluster variable.

---

## xml.R

### Objective

Generate the BEAST XML with a log-linear GLM substitution model for county-level phylogeographic inference.

### What it does

1. Similar to `DTACounty/scripts/xml.R` but replaces the standard substitution model with a `glmSubstitutionModel`.
2. The GLM models log transition rates as a linear combination of predictor variables:
   ```
   log(rate_ij) = beta_1 * X_1 + beta_2 * X_2 + ... + beta_n * X_n
   ```
   where each predictor has both an origin and destination version (e.g., sample size of origin county, sample size of destination county).
3. Includes a **design matrix** encoding predictor values for each county pair.
4. Uses **BSSVS** with coefficient indicator variables to perform Bayesian model selection — predictors can be included or excluded from the model.
5. MCMC operators include: scale, delta-exchange, random walk, multivariate normal, bit-flip.

### Predictors tested

| # | Predictor | Description |
|---|-----------|-------------|
| 1 | Sample size — origin | Number of sequenced genomes from origin county |
| 2 | Sample size — destination | Number of sequenced genomes from destination county |
| 3 | Median household income — origin | Census-derived economic indicator for origin county |
| 4 | Median household income — destination | Census-derived economic indicator for destination county |
| 5 | Population density — origin | People per square mile in origin county |
| 6 | Population density — destination | People per square mile in destination county |

### BEAST run parameters

| Parameter | Value |
|-----------|-------|
| Chain length | 100,000,000 |
| Log every | 10,000 states |

### Outputs

| File | Description |
|------|-------------|
| `GLMcounty.xml` | BEAST XML with GLM substitution model |

---

## 02_runBeast.sh

Submits 5 independent BEAST runs via SLURM (2 CPUs, 2 GB, 720-hour limit).

---

## 03_logCombine.sh

### What it does

Uses `logcombiner` with:
- **Burn-in:** 20,000,000 states
- **Resample:** every 200,000 states
- Combines into `GLMcounty/combined/combined_All_clades.log`

---

## 04_glmSummary.sh

### Objective

Summarize the GLM results: extract coefficient estimates, compute HPD intervals, posterior probabilities, and Bayes Factors for each predictor.

### What it does

Runs `GLM_Summary.R`, which:

1. Reads `GLMcounty/combined/combined_All_clades.log`.
2. For each predictor variable:
   - Extracts the coefficient time series and corresponding indicator time series.
   - Computes the **conditional effect size** (median of non-zero coefficient values) and 95% HPD interval.
   - Computes the **posterior probability (pp)** = mean of the indicator variable (proportion of MCMC samples where the predictor was included).
   - Computes the **Bayes Factor**:
     ```
     BF = [posterior / (1 - posterior)] / [prior / (1 - prior)]
     ```
     where the prior inclusion probability is `1 - exp(log(0.5) / n)` with `n` = number of predictor variables.
3. Classifies significance:
   - **Positive** (coral): 95% HPD entirely above 0
   - **Negative** (cyan): 95% HPD entirely below 0
   - **Non-Significant** (gray): HPD spans 0
4. Uses a posterior probability threshold of 0.20 (approximately BF > 3) to flag supported predictors.
5. Generates two figures.

### Interpreting the results

**`effect_size.pdf` — Conditional effect sizes:**
- **X-axis:** Conditional effect size (log scale, centered at 0)
- **Y-axis:** Predictor variables
- **Points:** Median conditional effect size. Color indicates significance direction (coral = positive, cyan = negative, gray = non-significant).
- **Horizontal error bars:** 95% HPD interval. If the interval crosses the black vertical line at 0, the effect is not statistically significant.
- **Positive effect sizes** mean the predictor increases viral diffusion rates between counties; negative means it decreases them.

**`posterior_prob.pdf` — Posterior inclusion probabilities:**
- **X-axis:** Posterior probability (0 to 1)
- **Y-axis:** Predictor variables
- **Bars:** Posterior probability that the predictor is included in the model. Color indicates whether the predictor exceeds the BF > 3 threshold (pp > 0.20).
- **Dashed blue line:** pp = 0.20 threshold (BF ~ 3)
- Predictors with bars extending past the threshold have statistical support for influencing transmission rates.

**Key findings from the paper:**
- **Destination county sample size** had the highest inclusion probability (pp = 1.0, decisive support), with a positive effect — more sequenced genomes in the destination county correlate with higher inferred transmission rates (likely reflecting better surveillance rather than causation).
- **Origin county sample size** (pp = 0.614) and **origin county population size** (pp = 0.387) also showed support with positive effects — larger, better-sampled origin counties export more virus.
- **Population density, median household income, and death counts** did not show sufficient statistical support as independent predictors of inter-county viral diffusion.

### Outputs

| File | Description |
|------|-------------|
| `effect_size.pdf` | Conditional effect sizes with 95% HPD intervals |
| `posterior_prob.pdf` | Posterior inclusion probabilities with BF threshold |

---

## R dependencies

| Package | Used by |
|---------|---------|
| `tidyverse`, `plotly` | pca.R |
| `dplyr`, `tidyr`, `readr`, `stringr` | GLM_Summary.R |
| `ggplot2`, `gridExtra` | GLM_Summary.R |
| `magrittr`, `lubridate`, `ape` | xml.R |
