library(ape)
library(dplyr)

# Directory with trees and CSV results
tree_dir <- "./tre"

# Find all CSV files ending in _BD_FFNN.csv
csv_files <- list.files(tree_dir, pattern = "_BD_FFNN.csv$", full.names = TRUE)

# Prepare summary dataframe
summary_df <- data.frame(Clade = character(),
                         Size = integer(),
                         R_naught = numeric(),
                         Infectious_period = numeric(),
                         stringsAsFactors = FALSE)

for (csv_file in csv_files) {
  # Extract clade name (e.g., Clade_17 from Clade_17_BD_FFNN.csv)
  clade_name <- sub("_BD_FFNN.csv$", "", basename(csv_file))
  tree_file <- file.path(tree_dir, paste0(clade_name, ".tre"))

  # Read the cleaned Newick tree
  if (file.exists(tree_file)) {
    tre <- read.tree(tree_file)
    clade_size <- length(tre$tip.label)
  } else {
    clade_size <- NA
  }

  # Read R_naught and Infectious_period
  csv_data <- read.csv(csv_file)
  R0 <- csv_data$R_naught[1]
  infectious <- csv_data$Infectious_period[1]

  # Append to summary
  summary_df <- rbind(summary_df, data.frame(Clade = clade_name,
                                             Size = clade_size,
                                             R_naught = R0,
                                             Infectious_period = infectious))
}

# Write summary table
write.csv(summary_df, file = file.path(tree_dir, "clade_summary.csv"), row.names = FALSE)
