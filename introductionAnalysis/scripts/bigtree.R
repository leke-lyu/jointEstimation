library(magrittr)
library(dplyr)
library(treeio)
library(lubridate)
library(ggtree)

pathToTree <- "tree137.nexus"
focalArea <- "Texas"

# read tree
tre <- pathToTree %>% treeio::read.beast()
tb <- tre %>% as_tibble()
root <- tb %>% filter(parent == node) %>% pull(node)

# find introduction events
importNode <- c()
for(j in tb$node){
  currentRow <- tb %>% filter(node == j)
  to <- (currentRow %>% pull(location))
  if( to == focalArea){
    from <- tb %>% filter(node == (currentRow %>% pull(parent))) %>% pull(location)
    if(from != focalArea){
      importNode <- c(importNode, j)
    }
  }
}

p <- ggtree(tre, size=1)

p + geom_point2(aes(subset=(node %in% importNode)), shape=23, size=1, fill="#ffde17") 

mrsd=(tb$num_date %>% max() %>% as.numeric() %>% lubridate::date_decimal() %>% as.Date())

