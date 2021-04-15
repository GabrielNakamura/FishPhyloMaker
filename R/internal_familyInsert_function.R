#' Internal add species to family node
#'
#' @param phy phylogenetic tree
#' @param spp_to_add_round2 species to add in family round of insertions
#'
#' @return a list with length two 
#'
insert_allbasefamily <- function(phy = phylo_order, spp_to_add_round2 = spp_to_add_round2, data_round3 = data_exRound3){
  count <- 0
  species_to_genus2 <- vector(mode = "list")
  add_all_family <- spp_to_add_round2 
  while(length(add_all_family) >= 1){
    count <- count + 1
    family_name <- data[match(add_all_family[1], 
                              data$s), "f"]
    node_family <- which(c(phylo_order$tip.label, 
                           phylo_order$node.label) == family_name)
    phylo_order <- phytools::bind.tip(tree = phylo_order, 
                                      tip.label = add_all_family[1], where = node_family, 
                                      position = 0)
    spp_data <- 1:nrow(data_exRound2)
    names(spp_data) <- data_exRound2$s
    insert_spp <- treedata_modif(phy = phylo_order, 
                                 data = spp_data, warnings = F)$nc$data_not_tree
    genus_in_tree <- sub("_.*", "", phylo_order$tip.label)[match(sub("_.*", 
                                                                     "", insert_spp), sub("_.*", "", phylo_order$tip.label))][!is.na(sub("_.*", 
                                                                                                                                         "", phylo_order$tip.label)[match(sub("_.*", 
                                                                                                                                                                              "", insert_spp), sub("_.*", "", phylo_order$tip.label))])]
    if (length(genus_in_tree) >= 1) {
      species_to_genus <- insert_spp[which(is.na(insert_spp[match(sub("_.*", 
                                                                      "", insert_spp), genus_in_tree)]) == FALSE)]
      species_to_genus2[[count]] <- species_to_genus
      for (i in 1:length(species_to_genus)) {
        phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                      species = species_to_genus[i])
      }
      insert_spp <- treedata_modif(phy = phylo_order, 
                                   data = spp_data, warnings = F)$nc$data_not_tree
    }
    rank_family2 <- unique(as.character(data[match(insert_spp, 
                                                   as.character(data$s)), 2]))
    if (length(rank_family2) == 0) {
      data_final <- 1:length(as.character(data$s))
      names(data_final) <- as.character(data$s)
      tree_res <- suppressWarnings(ape::drop.tip(phy = phylo_order, 
                                                 tip = treedata_modif(phy = phylo_order, 
                                                                      data = data_final)$nc$tree_not_data))
      data_exRound3 <- NULL
      res <- vector(mode = "list", length = 2)
      res[[1]] <- data_exRound3
      res[[2]] <- phylo_order
      return(res)
      break
    }
    else {
      add_all_family <- setdiff(insert_spp, 
                                   data_round3$s)
    }
  }
  res <- vector(mode = "list", length = 2)
  res[[1]] <- data_exRound3
  res[[2]] <- phylo_order
  return(res)
}
