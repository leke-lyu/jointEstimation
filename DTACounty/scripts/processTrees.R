library(dplyr)
library(magrittr)
library(treeio)
library(ggtree)

# Function to extract 'from' and 'to' locations
extractFromAndTo <- function(df){
  df <- df %>% mutate(to = location)
  df <- df %>% left_join(df %>% select(node, location), by = c("parent" = "node")) %>% 
    rename(from = location.y, location = location.x)
  return(df)
}

# Function to process a single tree file
process_tree <- function(tree_path) {
  tre <- read.beast(tree_path)
  tb <- fortify(tre)
  tb %<>% mutate(path = tree_path)
  tb %<>% select(parent, node, branch.length, location, isTip, path)
  tb %<>% extractFromAndTo()
  return(tb)
}

# Function to caculate lasting time
lasting_time <- function(n, tb, B = 0) {
  row <- tb %>% filter(node == n)
  p <- row %>% pull(parent)
  t <- row %>% pull(to)
  f <- row %>% pull(from)
  bl <- row %>% pull(branch.length)
  
  if (p == n) {
    return(B)
  } else {
    if (t == f) {
      B <- B + bl
      return(lasting_time(p, tb, B))
    } else {
      B <- B + bl / 2
      return(B)
    }
  }
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if there are at least two arguments (tree paths and output file path)
if (length(args) < 2) {
  stop("Please provide at least one tree file path and an output file path.")
}

# Get the output file path (last argument) and tree file paths (all preceding arguments)
output_file <- args[length(args)]
tree_paths <- args[-length(args)]

# Process each tree file and concatenate the results
all_fromAndTo <- tree_paths %>% 
  lapply(process_tree) %>% 
  bind_rows()


# Initialize the main dataframe
lt_df <- data.frame(Location = character(), LastingTime = numeric(), stringsAsFactors = FALSE)

# Loop through each unique path
for(i in unique(all_fromAndTo$path)) {
  tb <- all_fromAndTo %>% filter(path == i)
  Nodes <- tb %>% filter(isTip == TRUE) %>% pull(node)
  Locations <- tb %>% filter(isTip == TRUE) %>% pull(location)
  LastingTime <- c()
  
  # Calculate LastingTime for each node
  for(j in Nodes) {
    LastingTime <- c(LastingTime, lasting_time(j, tb))
  }
  
  # Create a temporary dataframe for this iteration
  temp_df <- data.frame(Location = Locations, LastingTime = LastingTime, stringsAsFactors = FALSE)
  
  # Concatenate the temporary dataframe with the main dataframe
  lt_df <- rbind(lt_df, temp_df)
}

df <- lt_df %>%
  group_by(Location) %>%
  summarise(MedianLastingTime = median(LastingTime, na.rm = TRUE))

# Source Sink Score
sss <- c()
districts <- (all_fromAndTo$to %>% unique() %>% sort())
for(i in districts){
  to <- all_fromAndTo %>%
    filter(from != i & to == i)
  from <- all_fromAndTo %>%
    filter(from == i & to != i)
  sss <- c(sss, (nrow(from) - nrow(to))/(nrow(from) + nrow(to))) 
}
names(sss) <- districts

df <- df %>%
  mutate(sss = sss[Location])

# Local Import Score
lis <- c()
for(i in districts){
  to <- all_fromAndTo %>%
    filter(from != i & to == i)
  mylocal <- all_fromAndTo %>%
    filter(from == i & to == i)
  lis <- c(lis, (nrow(to))/(nrow(mylocal) + nrow(to))) 
}
names(lis) <- districts

df <- df %>%
  mutate(lis = lis[Location])

# Save the concatenated results to an R data file
save(df, file = output_file)
