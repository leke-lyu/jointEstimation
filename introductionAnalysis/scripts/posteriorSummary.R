library(ggplot2)
library(lubridate)
library(tidyr)
library(magrittr)
library(dplyr)
library(ggpubr)
library(MASS)

directory <- "tree"
rdata_files <- list.files(path = directory, pattern = "\\.RData$", full.names = TRUE)
treeID <- rdata_files %>% basename(.) %>% sub("\\.RData$", "", .)

wholeList <- list()
for (file in rdata_files) {
  env <- new.env()
  load(file, envir = env)
  var_name <- ls(envir = env)
  data <- get(var_name, envir = env)
  wholeList <- append(wholeList, list(data))
}

wholeDF <- bind_rows(lapply(seq_along(wholeList), function(i) {
  wholeList[[i]] %>%
    mutate(whichTree = treeID[i])
}))

wholeDF %<>% filter(clusterSize != 0)

wholeDF <- wholeDF %>%
  mutate(category = case_when(
    clusterSize == 1 ~ "singletons",
    clusterSize > 1 & clusterSize <= 10 ~ "smaller than 10",
    clusterSize > 10 ~ "larger than 10"
  ))
wholeDF$date <- ((as.numeric(sub(".*-", "", wholeDF$epiweek))-1)/52 + as.numeric(sub("-.*", "", wholeDF$epiweek))) %>% date_decimal()
save(wholeDF, file = "wholeDF.RData")

count_df <- wholeDF %>%
  group_by(whichTree) %>%
  summarise(count = n())
print("number of introduction: ")
print(paste0("median: ", median(count_df$count)))
print(paste0("max: ", max(count_df$count)))
print(paste0("min: ", min(count_df$count)))
median_count <- median(count_df$count)
print("tree with median number of introduction: ")
count_df %>% filter(count == median(count_df$count)) %>% .$whichTree %>% print(.)

count_df <- wholeDF %>%
  group_by(whichTree, category) %>%
  summarise(count = n(), .groups = "drop") %>%
  filter(category == "singletons")
print("number of singletons: ")
print(paste0("median: ", median(count_df$count)))
print(paste0("max: ", max(count_df$count)))
print(paste0("min: ", min(count_df$count)))


count_df <- wholeDF %>%
  filter(source == "International") %>%
  group_by(whichTree) %>%
  summarise(count = n())
print("number of international introduction: ")
print(paste0("median: ", median(count_df$count)))
print(paste0("max: ", max(count_df$count)))
print(paste0("min: ", min(count_df$count)))

count_df <- wholeDF %>%
  filter(source == "Domestic") %>%
  group_by(whichTree) %>%
  summarise(count = n())
print("number of domestic introduction: ")
print(paste0("median: ", median(count_df$count)))
print(paste0("max: ", max(count_df$count)))
print(paste0("min: ", min(count_df$count)))

printIntroductionTimeBySize <- function(wholeDF){
  addDF <- wholeDF
  addDF$category <- "all"
  wholeDF <- bind_rows(wholeDF, addDF)
  
  weekly_counts <- wholeDF %>%
    group_by(date, whichTree, category) %>%
    summarise(introductions = n(), .groups = 'drop')
  
  all_trees <- unique(wholeDF$whichTree)
  all_date <- unique(wholeDF$date)
  all_category <- unique(wholeDF$category)
  
  weekly_counts <- expand.grid(whichTree = all_trees, date = all_date, category = all_category) %>%
    left_join(weekly_counts, by = c("whichTree", "date", "category")) %>%
    replace_na(list(introductions = 0))
  
  weekly_summary <- weekly_counts %>%
    group_by(date, category) %>%
    summarise(
      median = median(introductions),
      lower_bound = min(introductions),
      upper_bound = max(introductions),
      .groups = 'drop')
  
  weekly_summary <- weekly_summary %>%
    mutate(category = gsub("smaller than 10", "Clusters with size <= 10", category))
  
  weekly_summary <- weekly_summary %>%
    mutate(category = gsub("larger than 10", "Clusters with size > 10", category))
  
  weekly_summary <- weekly_summary %>%
    mutate(category = gsub("singletons", "Singletons", category))
  
  weekly_summary <- weekly_summary %>%
    mutate(category = gsub("all", "Total Imports", category))
  
  legend_order <- c("Total Imports", "Singletons", "Clusters with size <= 10", "Clusters with size > 10")
  colors <- c("Total Imports" = "grey", "Singletons" = "#038DB2", "Clusters with size <= 10" = "#FBB45C", "Clusters with size > 10" = "#81D0BB")
  
  p <- ggplot(weekly_summary, aes(x = date, color = category, fill = category)) +
    geom_line(aes(y = median)) +
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.4, color = NA) +
    scale_color_manual(values = colors, breaks = legend_order) +
    scale_fill_manual(values = colors, breaks = legend_order) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),  # Remove x-axis label
      axis.title.y = element_blank(),  # Remove y-axis label
      legend.position = "top",        # Move legend to the top
      legend.text = element_text(size = 5),      # Adjust legend text size
      legend.title = element_blank()     # Adjust legend title size
    )
  
  count <- weekly_summary %>% filter(date == unique(weekly_summary$date)[27]) %>% pull(median)
  paste0("In ", unique(weekly_summary$date)[27]) %>% print()
  paste0("size larger than 10: ", round((count[2]/count[1])*100, 1), " percent") %>% print()
  paste0("size smaller than 10: ", round((count[4]/count[1])*100, 1), " percent") %>% print()
  paste0("singletons: ", round((count[3]/count[1])*100, 1), " percent") %>% print()
  
  count <- weekly_summary %>% filter(date == unique(weekly_summary$date)[40]) %>% pull(median)
  paste0("In ", unique(weekly_summary$date)[40]) %>% print()
  paste0("size larger than 10: ", round((count[2]/count[1])*100, 1), " percent") %>% print()
  paste0("size smaller than 10: ", round((count[4]/count[1])*100, 1), " percent") %>% print()
  paste0("singletons: ", round((count[3]/count[1])*100, 1), " percent") %>% print()
  
  return(p)
}
p1 <- printIntroductionTimeBySize(wholeDF)

printIntroductionTimeBySource <- function(wholeDF){
  weekly_counts <- wholeDF %>%
    group_by(date, whichTree, source) %>%
    summarise(introductions = n(), .groups = 'drop')
  
  all_trees <- unique(wholeDF$whichTree)
  all_date <- unique(wholeDF$date)
  all_source <- unique(wholeDF$source)
  
  weekly_counts <- expand.grid(whichTree = all_trees, date = all_date, source = all_source) %>%
    left_join(weekly_counts, by = c("whichTree", "date", "source")) %>%
    replace_na(list(introductions = 0))
  
  weekly_summary <- weekly_counts %>%
    group_by(date, source) %>%
    summarise(
      median = median(introductions),
      lower_bound = min(introductions),
      upper_bound = max(introductions),
      .groups = 'drop')
  
  p <- ggplot(weekly_summary, aes(x = date, color = source, fill = source)) +
    geom_line(aes(y = median)) +
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha = 0.4, color = NA) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),  # Remove x-axis label
      axis.title.y = element_blank(),  # Remove y-axis label
      legend.position = "top",       # Move legend to the top
      legend.text = element_text(size = 6), # Adjust legend text size
      legend.title = element_blank()    # Adjust legend title size
    )
  return(p)
}
p2 <- printIntroductionTimeBySource(wholeDF)

ggsave("sizeSummary.pdf", p1, width = 8.9, height = 5, units = "cm")
ggsave("DvsI.pdf", p2, width = 8.9, height = 5, units = "cm")


parameters <- data.frame(
  whichTree = character(),
  k1 = numeric(),
  R1 = numeric(),
  k2 = numeric(),
  R2 = numeric(),
  k3 = numeric(),
  R3 = numeric(),
  stringsAsFactors = FALSE
)
for (tre in treeID) {
  # Filter data for the current tree
  treDF <- wholeDF %>% filter(whichTree == tre)
  
  # Filter Domestic and International sources
  doDF <- treDF %>% filter(source == "Domestic")
  inDF <- treDF %>% filter(source == "International")
  
  # Get cluster sizes
  size1 <- treDF$clusterSize
  size2 <- doDF$clusterSize
  size3 <- inDF$clusterSize
  
  # Initialize variables for k and R values
  k1 <- R1 <- k2 <- R2 <- k3 <- R3 <- NA
  
  # Fit negative binomial for treDF
  if (length(size1) > 0 && all(size1 >= 0)) {  # Ensure valid data for fitting
    fit1 <- fitdistr(size1, "negative binomial")
    k1 <- fit1$estimate["size"]
    R1 <- mean(size1)
  }
  
  # Fit negative binomial for doDF
  if (length(size2) > 0 && all(size2 >= 0)) {  # Ensure valid data for fitting
    fit2 <- fitdistr(size2, "negative binomial")
    k2 <- fit2$estimate["size"]
    R2 <- mean(size2)
  }
  
  # Fit negative binomial for inDF
  if (length(size3) > 0 && all(size3 >= 0)) {  # Ensure valid data for fitting
    fit3 <- fitdistr(size3, "negative binomial")
    k3 <- fit3$estimate["size"]
    R3 <- mean(size3)
  }
  
  # Append results to the results data frame
  parameters <- rbind(
    parameters,
    data.frame(
      whichTree = tre,
      k1 = k1, R1 = R1,
      k2 = k2, R2 = R2,
      k3 = k3, R3 = R3,
      stringsAsFactors = FALSE
    )
  )
}
summary_parameters <- parameters %>%
  summarise(
    k1_median = median(k1, na.rm = TRUE),
    k1_max = max(k1, na.rm = TRUE),
    k1_min = min(k1, na.rm = TRUE),
    
    R1_median = median(R1, na.rm = TRUE),
    R1_max = max(R1, na.rm = TRUE),
    R1_min = min(R1, na.rm = TRUE),
    
    k2_median = median(k2, na.rm = TRUE),
    k2_max = max(k2, na.rm = TRUE),
    k2_min = min(k2, na.rm = TRUE),
    
    R2_median = median(R2, na.rm = TRUE),
    R2_max = max(R2, na.rm = TRUE),
    R2_min = min(R2, na.rm = TRUE),
    
    k3_median = median(k3, na.rm = TRUE),
    k3_max = max(k3, na.rm = TRUE),
    k3_min = min(k3, na.rm = TRUE),
    
    R3_median = median(R3, na.rm = TRUE),
    R3_max = max(R3, na.rm = TRUE),
    R3_min = min(R3, na.rm = TRUE)
  )
print(summary_parameters)




