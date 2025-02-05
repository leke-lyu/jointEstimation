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
cleanTheTre <- function(myNode, myParent, tb, newtb, adjust){
  
  currentRow <- tb %>% filter(node == myNode)
  childrens <- tb %>% filter(parent == myNode) %>% pull(node)
  
  if(length(childrens) == 1){
    adjust <- adjust + (currentRow %>% pull(branch.length))
    newtb <- cleanTheTre(childrens, myParent, tb, newtb, adjust)
  }else{
    currentRow$parent <- myParent
    currentRow$branch.length <-  currentRow$branch.length + adjust
    newtb <- rbind(newtb, currentRow)
    if(length(childrens) == 2){
      for(child in childrens){
        newtb <- cleanTheTre(child, myNode, tb, newtb, 0)
      }
    }
  }
  return(newtb)
}

pathToTree <- "tree/tree137.nexus"
focalArea <- "Texas"
minimumSize <- 10

# read tree
tb <- pathToTree %>% treeio::read.beast() %>% as_tibble()
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

# select cluster by minimum size
keep_indices <- which(importSize > minimumSize)
importNode <- importNode[keep_indices]
importSize <- importSize[keep_indices]
treList <- treList[keep_indices]
tipList <- tipList[keep_indices]

# edit clusters
newtreList <- list()
for (i in seq_along(treList)) {
  tb <- treList[[i]]
  previous_nrow <- nrow(tb)
  repeat {
    tb <- tb %>% filter(!(is.na(label) & !(node %in% parent)))
    current_nrow <- nrow(tb)
    if (current_nrow == previous_nrow) break
    previous_nrow <- current_nrow
  }
  
  newtb <- tibble()
  newtb <- cleanTheTre((tb[1,] %>% pull(node)), (tb[1,] %>% pull(parent)), tb, newtb, 0)

  newtb <- newtb %>% mutate(node = if_else(!is.na(label), cumsum(!is.na(label)), as.integer(node)))
  
  count_non_na <- newtb %>% filter(!is.na(label)) %>% nrow()
  adjustment_index <- count_non_na + 1
  for (j in seq_len(nrow(newtb))) {
    if (is.na(newtb$label[j])) {
      current_node <- newtb$node[j]
      newtb$node[j] <- adjustment_index
      newtb$parent[newtb$parent == current_node] <- adjustment_index
      adjustment_index <- adjustment_index + 1
    }
  }
  
  rootEdge <- (newtb[1,] %>% pull(branch.length))
  newtb <- newtb[-1, ]
  t <- as.phylo(newtb)
  t$Nnode <- count_non_na - 1
  t$root.edge <- rootEdge
  newtreList[[i]] <- t
}

# test
writeTree <-function(tre, path){
  tree_content <- ape::write.tree(tre)
  nexus_content <- paste0(
    "#NEXUS\n",
    "begin trees;\n",
    '\ttree tree_1 = [&R] ', tree_content, "\n",
    "end;"
  )
  writeLines(nexus_content, path)
}
writeAllTrees <- function(treList, prefix = "t", extension = ".tre") {
  for (i in seq_along(treList)) {
    # Construct file name using the prefix and index
    file_name <- paste0(prefix, i, extension)
    
    # Write the tree to the file
    writeTree(treList[[i]], file_name)
  }
}
#writeAllTrees(newtreList)
cluster <- list(tipList, newtreList)
save(tipList, newtreList, file = "clusters.RData")



