#' Title Calculate the amount of phylogenetic deficit in assemblages
#'
#' @param phylo Phylogenetic tree in newick format, can be an object from [FishPhyloMaker()] function
#' @param data A data frame containing the classification  informing the level of insertions. This can be obtained
#'     from [FishPhyoMaker()] function
#' @param level Character indicating which level must be considered in the calculation of PD deficit.
#'     default is "Congeneric_insertion"
#'
#' @return A scalar containing the value of PD deficit for the level chosen
#' @export
#'
#' @seealso [FishPhyloMaker()] for phylogeny and data frame containing the classification of insertions
#' 
PD_defict <- function(phylo, data, level = "Congeneric_insertion"){
  names_exclude <- phylo$tip.label[na.omit(match(data[which(data$insertions == "Present_in_Tree"), "s"], 
                                                 phylo$tip.label))]
  
  if(all(!is.na(match( phylo$tip.label, names_exclude))) == TRUE){
    PD_present <- sum(phylo$edge.length)
    Darwinian_defict <- 0
  } else{
    exclude <- ape::drop.tip(phylo, tip = names_exclude)$tip.label
    phylo_present <- ape::drop.tip(phylo, tip = exclude)
    PD_present <- sum(phylo_present$edge.length)
    if(length(level) == 1){
      level_exclude <- ape::drop.tip(phylo, tip = phylo$tip.label[na.omit(match(data[which(data$insertions == level), "s"], 
                                                                                phylo$tip.label))])$tip.label
    }
    if(length(level) > 1){
      level_exclude <- ape::drop.tip(phylo, data$insertions[!is.na(match(data$insertions, level)), "s"])$tip.label 
    }
    phylo_level <- ape::drop.tip(phylo, level_exclude)
    PD_level <- sum(phylo_level$edge.length)
    PD_total <- sum(PD_present, PD_level)
    Darwinian_defict <- PD_level/PD_total
  }
  return(Darwinian_defict)
}
