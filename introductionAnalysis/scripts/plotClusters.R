library(ggtree)
library(dplyr)
library(ggtree)
library(patchwork)
library(ggplot2)

load("filteredHoustonDeltaMeta.RData")
load("clusters.RData")
custom_colors <- c(
  "Austin County" = "black",
  "Brazoria County" = "#FBB45C",
  "Chambers County" = "#164B6E",
  "Fort Bend County" = "#F9637C",
  "Galveston County" = "#FE7966",
  "Harris County" = "#038DB2",
  "Liberty County" = "#81D0BB",
  "Montgomery County" = "#8FA2DC",
  "Waller County" = "#FEAAC2"
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
  
    traits <- data.frame(GISAID_name = tre$tip.label) %>%
    left_join(filteredHoustonDeltaMeta, by = "GISAID_name") %>%
    pull(county)
  
  # Create the trait data frame
  trait_data <- data.frame(tip = tips, location = traits)
  
  # Create the ggtree plot
  p <- ggtree(tre, mrsd = t, size = 0.3) %<+% trait_data + 
    geom_tippoint(aes(color = location), shape = 21, size = 0.5, fill = "white", stroke = 0.8) +
    theme_tree2() +
    scale_x_continuous(name = "Time (Years)") +
    xlim(2021.2, 2021.8) +
    scale_color_manual("Location", values = custom_colors) +
    theme(
      legend.position = "none"
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

p1 <- create_tree_plot(whichTree = 78)
p2 <- create_tree_plot(whichTree = 53)
p3 <- create_tree_plot(whichTree = 45)
p4 <- create_tree_plot(whichTree = 54)
p5 <- create_tree_plot(whichTree = 10, TRUE)
plots <- list(p1, p2, p3, p4, p5)
combined_plot <- wrap_plots(plots, ncol = 1, heights = c(3, 2, 1.5, 0.5, 0.5))

ggsave("clusters.pdf", combined_plot, width = 6, height = 17.8, units = "cm")

dummy_data <- data.frame(
  location = names(custom_colors),
  x = 1:length(custom_colors),
  y = 1:length(custom_colors)
)

dummy_plot <- ggplot(dummy_data, aes(x = x, y = y, color = location)) +
  geom_point(shape = 21, size = 1, fill = "white", stroke = 1) +
  scale_color_manual("Location", values = custom_colors) +
  theme_void() +
  theme(legend.position = "bottom") # Place legend at the bottom


ggsave("dummy_legend.pdf", plot = dummy_plot, width = 20, height = 20, units = "cm")








