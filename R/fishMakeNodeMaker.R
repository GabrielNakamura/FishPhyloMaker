#' Title
#'
#' @param phylo A phylo object
#' @param data A three column data frame containing the name of species and their family and order. 
#'     can be obtained using 'FishTaxaMaker' function
#'
#' @return A phylo objetct with node names
#' @export
#'
FishNodeNameMaker <- function (phylo, data) 
{
  orders <- levels(as.factor(data$o))
  families <- levels(as.factor(data$f))
  data[match(orders[1], data$o), "s"]
  for (i in 1:length(orders)) {
    phylo <- ape::makeNodeLabel(phylo, "u", 
                                nodeList = list(Ord_name = data[match(orders[i], 
                                                                      data$o), "s"])
    )
    
    phylo$node.label[which(phylo$node.label == 
                             "Ord_name")] <- orders[i]
  }
  
  for (i in 1:length(families)) {
    phylo <- ape::makeNodeLabel(phylo, "u", 
                                nodeList = list(Ord_name = data[match(families[i], 
                                                                      data$f), "s"])
    )
    
    phylo$node.label[which(phylo$node.label == 
                             "Ord_name")] <- families[i]
  }
  
  return(phylo)
}
 


