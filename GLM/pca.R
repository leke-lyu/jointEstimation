# ============================================
# Title: PCA and Correlation-Based Variable Selection
# Purpose: Perform PCA to assess collinearity among epidemiological predictors;
#          then select 1 variable from a collinear cluster (case, death, population, samples)
#          based on lowest correlation with two fixed variables (MHI and PD).
# Author: [Leke Lyu]
# Date: [2025-06-04]
# ============================================

# -----------------------------
# Step 1: Load Required Libraries
# -----------------------------
library(tidyverse)
library(plotly)

# -----------------------------
# Step 2: Define Input Data
# -----------------------------
df <- tibble::tibble(
  county = c("Austin_County", "Brazoria_County", "Chambers_County", "Fort_Bend_County",
             "Galveston_County", "Harris_County", "Liberty_County", "Montgomery_County", "Waller_County"),
  case_count = c(2140, 36372, 4606, 72153, 40631, 398729, 7280, 50725, 3886),
  death_count = c(43, 625, 48, 603, 505, 6152, 265, 855, 66),
  population = c(30167, 372031, 46571, 822779, 350682, 4731145, 91628, 620443, 56794),
  sample_count = c(22, 192, 19, 568, 169, 5019, 45, 387, 34),
  median_household_income = c(68630, 86197, 102146, 105583, 75947, 68748, 61880, 95241, 75223),
  population_density = c(46.7, 274.0, 78.0, 955.1, 926.8, 2777.3, 79.1, 595.6, 110.6)
)

# -----------------------------
# Step 3: Perform PCA on All Predictors
# -----------------------------
df_numeric <- df %>% select(-county)
df_scaled <- scale(df_numeric)
pca <- prcomp(df_scaled, center = FALSE, scale. = FALSE)

# Extract loadings (PC1–PC3)
loadings <- as.data.frame(pca$rotation[, 1:3])
loadings$Variable <- rownames(loadings)

# -----------------------------
# Step 4: Visualize PCA Loadings (PC1 vs PC2 vs PC3)
# -----------------------------
plot_ly(loadings,
        x = ~PC1, y = ~PC2, z = ~PC3,
        text = ~Variable,
        type = "scatter3d", mode = "text+markers",
        marker = list(size = 5)) %>%
  layout(title = "PCA Loadings: PC1 vs PC2 vs PC3",
         scene = list(
           xaxis = list(title = "PC1"),
           yaxis = list(title = "PC2"),
           zaxis = list(title = "PC3")
         ))

# -----------------------------
# Step 5: Define Variable Groups
# -----------------------------
cluster_vars <- c("case_count", "death_count", "population", "sample_count")  # Variables from which we select one
fixed_vars <- c("median_household_income", "population_density")                                 # Variables to always retain

# -----------------------------
# Step 6: Compute Correlation Matrix
# -----------------------------
cor_mat_all <- cor(df_numeric[, c(cluster_vars, fixed_vars)], use = "complete.obs")

# -----------------------------
# Step 7: Define Scoring Function
# -----------------------------
# This function scores each candidate variable by its maximum absolute correlation
# with the fixed variables (MHI and PD). Lower scores = less collinear = preferred.
score_against_fixed <- function(cor_mat, candidates, fixed) {
  scores <- sapply(candidates, function(v) {
    max(abs(cor_mat[v, fixed]))
  })
  return(sort(scores))
}

# -----------------------------
# Step 8: Select Variable from Cluster
# -----------------------------
cor_scores_fixed <- score_against_fixed(cor_mat_all, candidates = cluster_vars, fixed = fixed_vars)
cat("Cluster variables ranked by max |r| with MHI and PD:\n")
print(cor_scores_fixed)

# Select the variable with lowest max correlation to fixed predictors
best_var <- names(cor_scores_fixed)[1]
cat("\nRecommended variable to keep from cluster:", best_var, "\n")

# -----------------------------
# Step 9: Define Final Predictor Set
# -----------------------------
final_vars <- c(best_var, fixed_vars)
cat("\nFinal selected predictor set:\n")
print(final_vars)
