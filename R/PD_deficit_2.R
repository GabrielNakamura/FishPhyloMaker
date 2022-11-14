
PD_deficit_2 <- 
function (phylo, data, level = "Congeneric_insertion") 
{
  if (is.null(phylo) == TRUE) {
    stop("/n A phylogenetic tree must be provided")
  }
  if (class(phylo) != "phylo") {
    stop("/n A phylo object must be provided")
  }
  names_exclude <- phylo$tip.label[na.omit(match(data[which(data$insertion == 
                                                              "present in tree"), "s"], phylo$tip.label))]
  if (all(!is.na(match(phylo$tip.label, names_exclude))) == 
      TRUE) {
    PD_present <- sum(phylo$edge.length)
    PD_level <- 0
    PD_total <- sum(phylo$edge.length)
    Darwinian_deficit <- PD_level/PD_total
    res <- c(PD_present, PD_level, PD_total, Darwinian_deficit)
    names(res) <- c("PDintree", "PDdeficit", "PDtotal", "Darwinian_deficit")
    return(res)
  }
  else {
    exclude <- ape::drop.tip(phylo, tip = names_exclude)$tip.label
    phylo_present <- ape::drop.tip(phylo, tip = exclude)
    PD_present <- sum(phylo_present$edge.length)
    if (length(level) == 1) {
      level_exclude <- ape::drop.tip(phylo, tip = phylo$tip.label[na.omit(match(data[which(data$insertions == 
                                                                                             level), "s"], phylo$tip.label))])$tip.label
    }
    if (length(level) > 1) {
      level_exclude <- ape::drop.tip(phylo, data[!is.na(match(data$insertion, 
                                                              level)), "s"])$tip.label
    }
    phylo_level <- ape::drop.tip(phylo, level_exclude)
    PD_level <- sum(phylo_level$edge.length)
    PD_total <- sum(phylo$edge.length)
    Darwinian_deficit <- PD_level/PD_total
    res <- c(PD_present, PD_level, PD_total, Darwinian_deficit)
    names(res) <- c("PDintree", "PDdeficit", "PDtotal", "Darwinian_deficit")
    return(res)
  }
}