# devtools::install_github("leke-lyu/TTAT")
library(magrittr)
library(dplyr)
library(TTAT)
myTTAT <- function(columnName, howMany){
  pb <- txtProgressBar(min = 0, max = howMany, style = 3)
  p <- numeric(0)  # Initialize p to store results
  tree_size <- numeric(0)
  
  for (i in 1:howMany) {
    
    tre <- newtreList[[i]]
    tip <- tre$tip.label
    trait_groups_V <- data.frame(GISAID_name = tre$tip.label) %>%
      left_join(filteredHoustonDeltaMeta, by = "GISAID_name") %>%
      pull(columnName)
    data <- data.frame(tipname = tip, trait = trait_groups_V)
    
    tips_to_drop <- data$tip[data$trait == "Unknown"]
    if(!is.null(tips_to_drop)){
      tre <- ape::drop.tip(tre, tips_to_drop)
      data %<>% filter(!(tipname %in% tips_to_drop))
      tip <- tre$tip.label
    }
    
    if((length(unique(data$trait)) > 1) & (nrow(data) > 10)){
      node <- setdiff(tre$edge[,1], tre$edge[,2])
      nTip <- length(tip)
      traits <- levels(factor(data[,2]))
      frequencyList <- frequencyOfTraits(tre, data, nTip, traits, node)
      distribution <- sapply(as.character(node), function(x) frequencyList[[x]])
      res <- aiNode(tre, data, nTip, traits, node)
      
      rep <- 1000
      nul <- numeric(rep)
      for(j in 1:rep){
        newTrait <- sample(rep(1:length(traits), distribution))
        newTrait <- traits[newTrait]
        d <- data.frame(tip, trait = newTrait)
        nul[j] <- aiNode(tre, d, nTip, traits, node)
      }
      P <- sum(round(nul, 8) <= round(res, 8)) / rep
      p <- c(p, P)
      
    } else {
      p <- c(p, NA)
    }
    
    tree_size <- c(tree_size, length(tre$tip.label))
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close the progress bar
  close(pb)
  
  return(data.frame(cluster_size = tree_size, p_value = p))
  
}

load("data/clusters.RData")
load("data/filteredHoustonDeltaMeta.RData")
filteredHoustonDeltaMeta$sex[is.na(filteredHoustonDeltaMeta$sex)] <- "Unknown"

age_groups <- myTTAT("new_age_groups", length(newtreList))
sex <- myTTAT("sex", length(newtreList))
save(age_groups, sex, file = "ttat.RData")

library(ggplot2)
age_groups$category <- "Age Groups"
sex$category <- "Sex"
age_groups <- age_groups %>% filter(!is.na(p_value))
age_groups %>% filter(`p_value` < 0.05)
sex <- sex %>% filter(!is.na(p_value))
sex %>% filter(`p_value` < 0.05)
combined_data <- bind_rows(age_groups, sex)
combined_data$category <- factor(combined_data$category, levels = c("Age Groups", "Sex"))
p <- ggplot(combined_data, aes(x = cluster_size, y = p_value, color = category, shape = category)) +
  geom_point() +
  scale_x_log10() +  # Log scale for x-axis
  labs(x = "Cluster Size (log scale)", y = "P-value", color = "Category", shape = "Category") +
  theme_minimal() +
  scale_color_manual(values = c("Age Groups" = "#038DB2", "Sex" = "#81D0BB")) +
  scale_shape_manual(values = c("Age Groups" = 19, "Sex" = 17)) +  # Manual shapes
  theme(legend.position = "top") +
  geom_hline(yintercept = 0.05, color = "red", linetype = "dashed", linewidth = 0.3)  +
  theme(
    plot.title = element_text(size = 8),   # Smaller title font size
    axis.title.x = element_text(size = 8), # Smaller x-axis title font size
    axis.title.y = element_text(size = 8),  # Smaller y-axis title font size
    legend.text = element_text(size = 8),   # Smaller legend text
    legend.title = element_text(size = 8),  # Smaller legend title
    legend.key.size = unit(0.5, "cm") 
  )

ggsave("ttat.pdf", p, width = 10, height = 5, units = "cm")

tre <- newtreList[[24]]
tip <- tre$tip.label
trait_groups_V <- data.frame(GISAID_name = tre$tip.label) %>%
  left_join(filteredHoustonDeltaMeta, by = "GISAID_name") %>%
  pull("new_age_groups")
data <- data.frame(tipname = tip, trait = trait_groups_V)

node <- setdiff(tre$edge[,1], tre$edge[,2])
nTip <- length(tip)
traits <- levels(factor(data[,2]))
frequencyList <- frequencyOfTraits(tre, data, nTip, traits, node)
distribution <- sapply(as.character(node), function(x) frequencyList[[x]])
res <- aiNode(tre, data, nTip, traits, node)
  
rep <- 1000
nul <- numeric(rep)
for(j in 1:rep){
  newTrait <- sample(rep(1:length(traits), distribution))
  newTrait <- traits[newTrait]
  d <- data.frame(tip, trait = newTrait)
  nul[j] <- aiNode(tre, d, nTip, traits, node)
}
P <- sum(round(nul, 8) <= round(res, 8)) / rep
p1 <- ggplot(data.frame(value = nul), aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill = "skyblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = res, color = "red", linetype = "dashed", linewidth = 1) +
  labs(
    title = paste0("p.value: ", P, "; AI: ", res),
    x = "Association Index",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8),   # Smaller title font size
    axis.title.x = element_text(size = 8), # Smaller x-axis title font size
    axis.title.y = element_text(size = 8)  # Smaller y-axis title font size
  )
ggsave("nulldistribution.pdf", p1, width = 5, height = 5, units = "cm")

library(ggtree)
custom_colors <- c(
  "Infants and Children" = "#FBB45C",
  "Middle-aged Adults" = "#F9637C",
  "Young Adults" = "#038DB2",
  "Seniors" = "#8FA2DC",
  "Teenagers" = "#FEAAC2",
  "Male" = "black",
  "Female" = "#999999"
)
create_tree_plot <- function(whichTree, show_x_axis = FALSE) {
  # Extract the specified tree and tips
  tre <- newtreList[[whichTree]]
  tips <- tre$tip.label
  # Extract metadata
  t <- filteredHoustonDeltaMeta %>% 
    filter(GISAID_name %in% tips) %>% 
    pull(collection_date) %>% 
    max()
  
  trait_groups_V <- data.frame(GISAID_name = tre$tip.label) %>%
    left_join(filteredHoustonDeltaMeta, by = "GISAID_name") %>%
    pull(new_age_groups)
  
  # Create the trait data frame
  trait_data <- data.frame(tip = tips, new_age_groups = trait_groups_V)
  
  # Create the ggtree plot
  p <- ggtree(tre, mrsd = t, size = 0.3) %<+% trait_data + 
    geom_tippoint(aes(color = new_age_groups), shape = 21, size = 0.5, fill = "white", stroke = 0.8) +
    theme_tree2() +
    scale_x_continuous(name = "Time (Years)") +
    xlim(2021.3, 2021.8) +
    scale_color_manual("Location", values = custom_colors) +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 1),  # Adjust text size
      legend.title = element_text(size = 1), # Adjust title size
      legend.key.size = unit(0.5, "cm")       # Adjust key size
    )
  
  # Adjust x-axis display based on the input
  if (!show_x_axis) {
    p <- p + theme(
      axis.title.x = element_blank(),   # Remove x-axis title
      axis.text.x = element_blank(),    # Remove x-axis text
      axis.ticks.x = element_blank(),   # Remove x-axis ticks
      axis.line.x = element_blank()     # Remove x-axis line
    )
  }
  
  return(p)
}
p2 <- create_tree_plot(whichTree = 24, TRUE)
ggsave("mycluster.pdf", p2, width = 5, height = 5, units = "cm")

