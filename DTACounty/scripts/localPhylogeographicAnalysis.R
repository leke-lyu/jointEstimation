library(ggplot2)
library(dplyr)
library(stringr)
library(magrittr)
library(sf)
library(tigris)
library(viridis)

combine_rdata_files <- function(path) {
  # Ensure the path ends with a slash
  if (!grepl("/$", path)) path <- paste0(path, "/")
  
  # List all .RData files in the directory
  rdata_files <- list.files(path, pattern = "\\.RData$", full.names = TRUE)
  
  # Initialize an empty list to store data frames
  data_frames <- list()
  
  # Load and extract each .RData file
  for (file in rdata_files) {
    # Load the RData file
    env <- new.env()  # Use a separate environment to avoid overwriting variables
    load(file, envir = env)
    
    # Assuming the RData contains a single data frame
    data_frame <- get(ls(env)[1], envir = env)  # Get the first object in the environment
    
    # Check if the object is a data frame
    if (is.data.frame(data_frame)) {
      data_frames[[length(data_frames) + 1]] <- data_frame
    } else {
      warning(paste("Skipping", file, "as it does not contain a data frame"))
    }
  }
  
  # Combine all data frames using rbind
  combined_df <- do.call(rbind, data_frames)
  
  return(combined_df)
}
process_string <- function(s) {
  
  trimmed_string <- str_remove(s, "^location\\.rates\\.")
  
  match <- str_locate(trimmed_string, "_District\\.|_County\\.")
  
  if (!is.na(match[1])) {
    first_part <- str_sub(trimmed_string, 1, match[2] - 1) # Extract first part
    second_part <- str_sub(trimmed_string, match[2] + 1, str_length(trimmed_string)) # Extract second part
    return(list(first_part, second_part))
  } else {
    return(NULL)
  }
  
}
find_k <- function(m) {
  for (n in 1:m) {
    if (m == n * (n - 1)) {
      return(n)
    }
  }
  return(NULL)
}
getBF <- function(myLog){
  selected_Log <- myLog %>%
    select(starts_with("location.rates"), starts_with("location.indicators")) %>% 
    mutate(across(starts_with("location.rates"), ~ . * get(gsub("rates", "indicators", cur_column())), .names = "real.rates.{gsub('location.rates.', '', col)}"))
  k <- (length(colnames(selected_Log))/3) %>% find_k()
  colNames <- colnames(selected_Log)[1:(dim(selected_Log)[2]/3)]
  from <- lapply(colNames, process_string) %>% sapply(., function(x) x[[1]])
  to <- lapply(colNames, process_string) %>% sapply(., function(x) x[[2]])
  posterior.mean.of.indicators <- numeric(k*(k-1))
  for(i in 1:(k*(k-1))){
    posterior.mean.of.indicators[i] <- paste0("location.indicators.", from[i], ".", to[i]) %>% pull(selected_Log, .) %>% mean()
  }
  posterior.mean.of.location.rates <- numeric(k*(k-1))
  for(i in 1:(k*(k-1))){
    posterior.mean.of.location.rates[i] <- paste0("location.rates.", from[i], ".", to[i]) %>% pull(selected_Log, .) %>% mean()
  }
  posterior.mean.of.real.rates <- numeric(k*(k-1))
  for(i in 1:(k*(k-1))){
    posterior.mean.of.real.rates[i] <- paste0("real.rates.", from[i], ".", to[i]) %>% pull(selected_Log, .) %>% .[. != 0] %>% { if(length(.) == 0) 0 else mean(.) }
  }
  bf <- numeric(k*(k-1))
  q <- (log(2)+k-1)/(k*(k-1))
  for(i in 1:(k*(k-1))){
    p <- posterior.mean.of.indicators[i]
    bf[i] <- (p*(1-q))/((1-p)*q)
  }
  bf_df <- data.frame(from, to, posterior.mean.of.indicators, posterior.mean.of.location.rates, posterior.mean.of.real.rates, bf)
  
  return(bf_df)
}

path <- "output/combined"

result <- path %>% combine_rdata_files() %>%
  group_by(Location) %>%  # Group data by counties
  summarise(
    MedianLastingTime_min = min(MedianLastingTime, na.rm = TRUE),
    MedianLastingTime_max = max(MedianLastingTime, na.rm = TRUE),
    MedianLastingTime_median = median(MedianLastingTime, na.rm = TRUE),
    sss_min = min(sss, na.rm = TRUE),
    sss_max = max(sss, na.rm = TRUE),
    sss_median = median(sss, na.rm = TRUE),
    lis_min = min(lis, na.rm = TRUE),
    lis_max = max(lis, na.rm = TRUE),
    lis_median = median(lis, na.rm = TRUE)
  )
result$Location <- gsub("_County", "", result$Location)
result$Location <- gsub("_", " ", result$Location)
result <- result %>%
  mutate(Location = factor(Location, levels = c(
    "Harris", "Fort Bend", "Montgomery",
    "Brazoria", "Galveston", "Liberty",
    "Waller", "Chambers", "Austin"
  )))

selected_counties <- counties(state = "TX", cb = TRUE) %>% filter(NAME %in% result$Location)
p1 <- ggplot() +
  geom_sf(data = selected_counties, fill = "dodgerblue3", color = "white", linewidth = 1) +
  theme_minimal() +
  theme(legend.position = "none")

locations_data <- data.frame(
  regions = c("Harris", "Fort Bend", "Montgomery", "Brazoria", "Galveston", "Liberty", "Waller", "Chambers", "Austin"),
  latitude = c(29.699528, 29.534355, 30.315601, 29.17403, 29.39084, 30.12546, 29.98955, 29.7626, 29.91717),
  longitude = c(-95.418809, -95.750403, -95.529414, -95.48846, -95.03648, -94.89097, -95.98678, -94.62657, -96.32391))

pathToLog <- paste0(path, "/combined_All_clades.log")
burnin <- 0
myLog <- read.csv(pathToLog, sep = "\t", header = TRUE)
myLog <- myLog[(floor(burnin*0.01*(dim(myLog)[1]))+1):(dim(myLog)[1]),]
bf_df <- myLog %>% getBF()
bf_df %<>% filter(posterior.mean.of.indicators > 0.5)
bf_df %<>% rename(TransitionRate = posterior.mean.of.real.rates)
bf_df$from <- gsub("_County", "", bf_df$from)
bf_df$from <- gsub("_", " ", bf_df$from)
bf_df$to <- gsub("_County", "", bf_df$to)
bf_df$to <- gsub("_", " ", bf_df$to)
write.csv(bf_df, file = "output_tb.csv", row.names = FALSE)

bf_df <- merge(bf_df, locations_data, by.x = "from", by.y = "regions", all.x = TRUE)
names(bf_df)[names(bf_df) == "latitude"] <- "from_latitude"
names(bf_df)[names(bf_df) == "longitude"] <- "from_longitude"
bf_df <- merge(bf_df, locations_data, by.x = "to", by.y = "regions", all.x = TRUE)
names(bf_df)[names(bf_df) == "latitude"] <- "to_latitude"
names(bf_df)[names(bf_df) == "longitude"] <- "to_longitude"
bf_df %<>% filter( bf > 100)

p1 <- p1 + 
  geom_curve(data = bf_df, 
             aes(x = from_longitude, y = from_latitude, xend = to_longitude, yend = to_latitude, linewidth = TransitionRate),
             arrow = arrow(angle = 15, length = unit(0.2, "cm")),
             curvature = 0.1,
             color = "dimgray") +
  labs(x = NULL, y = NULL) +
  scale_linewidth_continuous(range = c(0.2, 2.0)) +
  theme_minimal() +
  theme(legend.position = "top")

p2 <- ggplot(result, aes(x = Location)) +
  geom_point(aes(y = sss_median), shape = 23, fill = "#2ca02c", size = 3) + # Blue diamond for median
  geom_errorbar(aes(ymin = sss_min, ymax = sss_max), color = "goldenrod", size = 0.3, width = 0.2) + # Golden error bars
  scale_y_continuous(limits = c(-1, 1)) + # Set y-axis range
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank()
  )

p3 <- ggplot(result, aes(x = Location)) +
  geom_point(aes(y = lis_median), shape = 22, fill = "#ff7f0e", size = 3) + # Red square for median
  geom_errorbar(aes(ymin = lis_min, ymax = lis_max), color = "goldenrod", size = 0.3, width = 0.2) + # Golden error bars
  scale_y_continuous(limits = c(0, 1)) + # Set y-axis range
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank()
  ) 

p4 <- ggplot(result, aes(x = Location)) +
  geom_point(aes(y = MedianLastingTime_median), shape = 21, fill = "#1f77b4", size = 3) + # Green circle for median
  geom_errorbar(aes(ymin = MedianLastingTime_min, ymax = MedianLastingTime_max), color = "goldenrod", size = 0.3, width = 0.2) + # Golden error bars
  scale_y_continuous(limits = c(0, 0.6)) + # Set y-axis range
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )


library(ggpubr)

metric <- ggarrange(p2, p3, p4,
                           nrow = 3, ncol = 1,
                           heights = c(1.25, 1.25, 2),
                           align = "v"
)
ggsave("metric.pdf", metric, width = 10, height = 10, units = "cm")
ggsave("map.pdf", p1, width = 9, height = 9, units = "cm")




