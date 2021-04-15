#' Internal add species to order node
#'
#' @param phylo phylogenetic tree
#' @param data_round3 data frame
#'
#' @return phylogenetic tree
#'
insert_allbaseOrd <- function(phylo = phylo_order, data_round3 = data_exRound3, list_names_ord = list_order){
  order_add_round3 <- unique(data_round3$o)
  for(i in 1:nrow(data_round3)){
    node_order_pos <- which(c(phylo$tip.label, 
                              phylo$node.label) == data_round3$o[i])
    phylo <- phytools::bind.tip(tree = phylo, 
                                      tip.label = data_round3$s[i], where = node_order_pos, 
                                      position = 0)
  }
  return(phylo)
}
