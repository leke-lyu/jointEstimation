library(magrittr)
library(dplyr)
library(treeio)
library(lubridate)
extractSubTre <- function(myNode, root, tb, subtreTB){
  
  currentRow <- tb %>% filter(node == myNode)
  loc <- currentRow %>% pull(location)
  
  if(loc == focalArea){
    subtreTB <- rbind(subtreTB, currentRow)
    childrens <- tb %>% filter(parent == myNode) %>% pull(node)
    if(myNode == root){
      childrens <- childrens[!childrens %in% root]
    }
    for(child in childrens){
      subtreTB <- extractSubTre(child, root, tb, subtreTB)
    }
  }
  
  return(subtreTB)
}
addHeight <- function(tb, node, height, root){
  childrens <- tb %>% filter(parent == !!node) %>% pull(node)
  if(node == root){
    childrens <- childrens[!childrens %in% root]
  }
  for (currentNode in childrens) {
    currentRow <- tb %>% filter(node == currentNode)
    height[currentNode] <- currentRow$branch.length + height[currentRow$parent]
    if((currentRow$label %>% is.na())){
      height <- addHeight(tb, currentNode, height, root)
    }
  }
  return(height)
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) stop("3 arguments must be provided: tree file and focal area and maxDate")
pathToTree <- args[1]
focalArea <- args[2]
maxDate <- args[3]
dir <- dirname(pathToTree)
treeName <- sub("\\.nexus$", "", basename(pathToTree))

# read tree
tb <- pathToTree %>% treeio::read.beast() %>% as_tibble()
root <- tb %>% filter(parent == node) %>% pull(node)

# add height as a vector
Height <- numeric(nrow(tb))
Height <- addHeight(tb, root, Height, root)
maxHeight <- max(Height)

# find introduction events
importNode <- c()
importSource <- c()
importHeight <- c()
for(j in tb$node){
  currentRow <- tb %>% filter(node == j)
  to <- (currentRow %>% pull(location))
  if( to == focalArea){
    from <- tb %>% filter(node == (currentRow %>% pull(parent))) %>% pull(location)
    if(from != focalArea){
      importNode <- c(importNode, j)
      importSource <- c(importSource, from)
      importHeight <- c(importHeight, (Height[j] - (currentRow %>% pull(branch.length))/2))
    }
  }
}

year <- (importHeight - maxHeight +  as.numeric(maxDate)) %>% date_decimal(.) %>% epiyear(.)
week <- (importHeight - maxHeight +  as.numeric(maxDate)) %>% date_decimal(.) %>% epiweek(.)
importEpiweek <- paste0(year, "-", week)

# extract clusters and get cluster size
treList <- list()
for(myNode in importNode){
  subtreTB <- tibble()
  subtreTB <- extractSubTre(myNode, root, tb, subtreTB)
  treList[[length(treList) + 1]] <- subtreTB
}
tipList <- lapply(treList, function(tb) { #tipVector <- unlist(tipList)
  tb %>% filter(!is.na(label)) %>% pull(label)
})
importSize <- sapply(tipList, length)

importSummary <- data.frame(node = importNode, source = importSource, epiweek = importEpiweek, clusterSize = importSize) 
assign(treeName, importSummary)
save(list=treeName, file = paste0(dir, "/", treeName, ".RData"))
