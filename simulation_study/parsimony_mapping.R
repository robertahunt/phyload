tree <- rtree(10,FALSE)

chars <- sample(1:4,10,replace=TRUE)

# X MUXT BE ORDERED IN ORDER THAT TIPS ARE IN TIP LABELS
siteFitchMappingUnrooted <- function(phy,X) {
  recover()

  ntax <- length(phy$tip.label)

  phy_root_children <- phy$edge[phy$edge[,1] == ntax+1,2]

  phy <- reorder(phy,"postorder")

  rooted_edge <- phy$edge

  new_root <- ntax + phy$Nnode + 1

  rooted_edge <- rbind(rooted_edge,c(new_root,ntax+1))

  rooted_edge[rooted_edge[,2] == phy_root_children[3],1] <- new_root

  states <- vector("list",phy$Nnode + 1)

  for (i in 1:length(X)) {
    states[[i]] <- X[i]
  }

  for (i in 1:(dim(rooted_edge)[1]/2)) {
    row2 <- 2*i
    row1 <- row2 - 1
    int <- intersect(states[[rooted_edge[row1,2]]],states[[rooted_edge[row2,2]]])
    if (length(int) > 0) {
      parent_states <- int
    } else {
      parent_states <- union(states[[rooted_edge[row1,2]]],states[[rooted_edge[row2,2]]])
    }
    states[[rooted_edge[row1,1]]] <- parent_states
    if () {

    }
  }

}



siteFitchMapping(tree,chars)
