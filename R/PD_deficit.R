#' Title Calculate the amount of phylogenetic deficit in assemblages
#'
#' @param phylo Phylogenetic tree in newick format, can be an object from \code{\link{FishPhyloMaker}} function
#' @param data A data frame containing the classification  informing the level of insertions. This can be obtained
#'     from \code{\link{FishPhyloMaker}} function
#' @param level Character indicating which level must be considered in the calculation of PD deficit. 
#'     Can be a vector with the levels ("Congeneric_insertion", "Congeneric_Family_level", "Family_insertion", "Order_insertion") 
#'     which will be considered in the calculation of phylogenetic deficit.
#'     default is "Congeneric_insertion".
#'
#' @return A vector containing four values: 
#' 
#'     - Amount phylogenetic information present in the tree before insertions (PDintree)
#'     
#'     - Amount of phylogenetic information inserted in the tree (PDdeficit)
#'     
#'     - Total Phylogenetic information of the tree (PDtotal)
#'     
#'     - A ratio calculated as PDdeficit/PDtotal (Darwinian_deficit)
#' @export
#'
#' @seealso \code{\link{FishPhyloMaker}} for phylogeny and data frame containing the classification of insertions
#' 
PD_deficit <- function(phylo, data, level = "Congeneric_insertion"){
  if(is.null(phylo) == TRUE){
    stop("/n A phylogenetic tree must be provided")
  }
  if(class(phylo) != "phylo"){
    stop("/n A phylo object must be provided")
  }
  names_exclude <- phylo$tip.label[na.omit(match(data[which(data$insertions == "Present_in_Tree"), "s"], 
                                                 phylo$tip.label))]
  
  if(all(!is.na(match(phylo$tip.label, names_exclude))) == TRUE){
    PD_present <- sum(phylo$edge.length)
    PD_level <- 0
    PD_total <- sum(PD_present, PD_level)
    Darwinian_deficit <- PD_level/PD_total
    res <- c(PD_present, PD_level, PD_total, Darwinian_deficit)
    names(res) <- c("PDintree", "PDdeficit", "PDtotal", "Darwinian_deficit")
    return(res)
  } else{
    exclude <- ape::drop.tip(phylo, tip = names_exclude)$tip.label
    phylo_present <- ape::drop.tip(phylo, tip = exclude)
    PD_present <- sum(phylo_present$edge.length)
    if(length(level) == 1){
      level_exclude <- ape::drop.tip(phylo, tip = phylo$tip.label[na.omit(match(data[which(data$insertions == level), "s"], 
                                                                                phylo$tip.label))])$tip.label
    }
    if(length(level) > 1){
      level_exclude <- ape::drop.tip(phylo, data[!is.na(match(data$insertions, level)), "s"])$tip.label 
    }
    phylo_level <- ape::drop.tip(phylo, level_exclude)
    PD_level <- sum(phylo_level$edge.length)
    PD_total <- sum(PD_present, PD_level)
    Darwinian_deficit <- PD_level/PD_total
    res <- c(PD_present, PD_level, PD_total, Darwinian_deficit)
    names(res) <- c("PDintree", "PDdeficit", "PDtotal", "Darwinian_deficit")
    return(res)
  }
}