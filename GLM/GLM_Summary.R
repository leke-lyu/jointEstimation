library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(dplyr)
library(gridExtra)

# Function to calculate the HPD interval for a given vector
calc_min_interval <- function(x, alpha = 0.05) {
  x <- sort(x)
  n <- length(x)
  cred_mass <- 1.0 - alpha
  
  interval_idx_inc <- floor(cred_mass * n)
  n_intervals <- n - interval_idx_inc
  interval_width <- x[(interval_idx_inc + 1):n] - x[1:n_intervals]
  
  if (length(interval_width) == 0) {
    stop("Too few elements for interval calculation")
  }
  
  min_idx <- which.min(interval_width)
  hdi_min <- x[min_idx]
  hdi_max <- x[min_idx + interval_idx_inc]
  
  return(c(hdi_min, hdi_max))
}

# Function to calculate the Bayes Factor for indicators
BayesFactor <- function(indicators, n) {
  prior <- 1 - exp(log(0.5) / n)
  posterior <- mean(indicators)
  if (posterior == 1) {
    posterior <- (length(indicators) - 1) / length(indicators)
  }
  BF <- (posterior / (1 - posterior)) / (prior / (1 - prior))
  return(BF)
}

pathToLog <- paste0("GLMcounty/combined/combined_All_clades.log")
burnin <- 0
trait <- "location"
varNum <- 10

myLog <- read.csv(pathToLog, sep = "\t", header = TRUE, comment.char = "#")
myLog <- myLog[(floor(burnin*0.01*(dim(myLog)[1]))+1):(dim(myLog)[1]),]

output <- data.frame(
  Variable = character(),
  CoefMedian = numeric(),
  CoefLowerHPD = numeric(),
  CoefUpperHPD = numeric(),
  pp = numeric(),
  BF = numeric(),
  stringsAsFactors = FALSE
)

for (col in names(myLog)) {
  if (grepl("Times", col)) {
    varID <- str_extract(col, "\\d+")
    indcol <- paste0(trait, ".coefIndicators", varID)
    indicators <- myLog[[indcol]]
    
    if (max(indicators) == 0) {
      interval <- c(0, 0)
      median <- 0
    } else {
      NoZero <- myLog[[col]]
      NoZero[NoZero == 0] <- NA
      NoNan <- na.omit(NoZero)
      interval <- calc_min_interval(NoNan)
      median <- median(NoNan, na.rm = TRUE)
    }
    
    pp <- mean(indicators)
    bf <- BayesFactor(indicators, as.numeric(varNum))
    
    output <- rbind(output, data.frame(
      Variable = varID,
      CoefMedian = median, CoefLowerHPD = interval[1], 
      CoefUpperHPD = interval[2], pp = pp, BF = bf
    ))
  }
}

# Prepare data for plotting
output <- output %>%
  mutate(
    Significant = case_when(
      CoefLowerHPD > 0 ~ "Positive",             # Statistically supported positive association (red)
      CoefUpperHPD < 0 ~ "Negative",             # Statistically supported negative association (blue)
      TRUE ~ "Non-Significant"                   # No association (gray)
    ),
    AboveThreshold = case_when(
      pp > 0.20 & Significant == "Negative" ~ "AboveNegative",  # For bars above threshold with negative effect
      pp > 0.20 ~ "Above",                                   # For bars above threshold without negative effect
      TRUE ~ "Below"                                         # Below threshold
    )
  )
output$Variable <- c("sample size - origin", "sample size - destination", "median household income - origin", "median household income - destination", "population density - origin", "population density - destination")

# Define custom colors based on significance and threshold
custom_colors <- c("Positive" = "coral", "Negative" = "cyan", "Non-Significant" = "gray")
threshold_colors <- c("AboveNegative" = "cyan", "Above" = "coral", "Below" = "lightgray")

# Part 1: Plot for Conditional Effect Size with HPD Intervals
effect_size_plot <- ggplot(output, aes(x = CoefMedian, y = Variable)) +
  geom_point(aes(color = Significant, fill = Significant), shape = 21, size =4) +  # Shape 21 allows filling
  geom_errorbarh(aes(xmin = CoefLowerHPD, xmax = CoefUpperHPD), height = 0.2, color = "black") +
  geom_vline(xintercept = 0, linetype = "solid", color = "black") +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  scale_size_continuous(range = c(1, 5), guide = "none") +
  scale_y_discrete(position = "right") +
  labs(
    x = "Conditional Effect Size", 
    y = "Predictor Variable",
    title = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 8),        # Make legend text smaller
    legend.key.size = unit(0.5, "lines"),       # Reduce size of legend keys
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 6, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 6)
  )

# Part 2: Plot for Posterior Probability with additional BF = 3 and BF = 100 threshold lines
posterior_prob_plot <- ggplot(output, aes(x = pp, y = Variable, fill = AboveThreshold)) +
  geom_col(color = "black", width = 0.6) +
  geom_vline(xintercept = 0.20, linetype = "dashed", color = "darkblue") +  # BF = 3 threshold line
  #geom_vline(xintercept = 1.0, linetype = "dashed", color = "darkblue") +  # BF = 100 threshold line (adjust as necessary)
  scale_fill_manual(values = threshold_colors) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_discrete(position = "right") +
  labs(
    x = "Posterior Probability",
    y = NULL,
    title = NULL
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 8),        # Make legend text smaller
    legend.key.size = unit(0.5, "lines"),       # Reduce size of legend keys
    axis.text.y = element_text(size = 6, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 6)
  )

# Arrange plots side-by-side

ggsave("posterior_prob.pdf", posterior_prob_plot, width = 7.8, height = 10, units = "cm")
ggsave("effect_size.pdf", effect_size_plot, width = 7.8, height = 10, units = "cm")
