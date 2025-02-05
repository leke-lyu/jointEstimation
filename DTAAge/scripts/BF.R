library(ggplot2)
library(dplyr)
library(stringr)
library(magrittr)
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
  trimmed_string <- str_remove(colNames, "^location\\.rates\\.")
  split_strings <- str_split(trimmed_string, "\\.", simplify = TRUE)
  from <- split_strings[, 1]
  to <- split_strings[, 2]
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
get_color <- function(x) {
  trait <- sub("^(src-|snk-)", "", x)
  traits_colors[trait]
}

pathToLog <- "output/combined/combined_All_clades.log"
burnin <- 0
myLog <- read.csv(pathToLog, sep = "\t", header = TRUE)
myLog <- myLog[(floor(burnin*0.01*(dim(myLog)[1]))+1):(dim(myLog)[1]),]
colnames(myLog) <- gsub("Middle\\.aged_Adults", "Middle-aged_Adults", colnames(myLog))

bf_df <- myLog %>% getBF()
bf_df %<>% filter(posterior.mean.of.indicators > 0.5)
bf_df %<>% rename(TransitionRate = posterior.mean.of.real.rates)
bf_df %<>%
  mutate(group = case_when(
    bf > 100 ~ "type_a",
    bf > 30 ~ "type_b",
    bf > 10 ~ "type_c",
    bf > 3 ~ "type_d",
    TRUE ~ "else"
  ))
  
save(bf_df, file = "bf.RData")

traits <- c(bf_df$from, bf_df$to) %>% unique()
traits_colors <- c("#038DB2", "#F9637C", "#FBB45C", "#81D0BB", "#8FA2DC")
names(traits_colors) <- traits
custom_colors <- c("type_a" = "grey10", "type_b" = "grey30", "type_c" = "grey60", "type_d" = "grey90")

library(ggplot2)
library(circlize)
df <- bf_df
df$from <- paste0("src-", df$from)
df$to <- paste0("snk-", df$to)
new_traits <- c("src-Teenagers", "src-Young_Adults", "src-Seniors", "src-Middle-aged_Adults", "src-Infants_and_Children", "snk-Seniors", "snk-Middle-aged_Adults", "snk-Teenagers", "snk-Infants_and_Children", "snk-Young_Adults")
grid.col <- setNames(sapply(new_traits, get_color), new_traits)
df <- df %>%
  mutate(bf_color = custom_colors[group])

all_sectors <- unique(c(df$from, df$to))
gap.after <- rep(1, length(all_sectors))
gap.after[length(unique(df$from))] <- 20
gap.after[length(all_sectors)] <- 20

circos.clear()
circos.par(start.degree = 270, gap.after = gap.after)

pdf("flow.pdf", width=7.8 / 2.54, height=7.8 / 2.54)
chordDiagram(df[,c("from", "to", "TransitionRate")],
             order=new_traits,
             grid.col = grid.col,
             col = df$bf_color,
             directional = 1,
             direction.type = c("arrows"),
             link.arr.type = "big.arrow",
             annotationTrack = "grid", 
             annotationTrackHeight = mm_h(5))

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "outside", niceFacing = TRUE, adj = c(0.5, 0.5), cex = 0.9, family = "serif")
}, bg.border = NA)

# Add a legend
legend("topright", legend = names(custom_colors), fill = custom_colors, cex = 0.5)

# Close the PDF
dev.off()





