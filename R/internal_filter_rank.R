## function to downloading phylogeny from all orders in data
filter_rank<- function(order){
  if(length(which(sub("_.*", "", unlist(order)) == "not.found")) >= 1){
    
    phy_ord<- fishtree::fishtree_phylogeny(unlist(order)[-which(sub("_.*", "", unlist(order)) == "not.found")])
    
  } else{
    phy_ord<- fishtree::fishtree_phylogeny(unlist(order))
  }
  phy_ord
}

