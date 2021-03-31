#' Function to inform which species must be added to the mega-tree phylogeny in the insertion process.
#'
#' @param data A data frame with three column containing the name of species (s), the Family (f) and Order (o). This 
#'     can be generated with function \code{\link{FishTaxaMaker}}
#'     
#' @details This function can be used  in order to known which species that must be added in the insertion process 
#'     made by \code{\link{FishPhyloMaker}}.
#'
#' @return A data frame containing a column informing at which level the species in data must be added.
#' 
#' @export
#'
#' @examples
#' \dontrun{
#'     data(neotropical_comm)
#'     data_comm <- neotropical_comm[, -c(1, 2)]
#'     taxon_data <- tab_function(data_comm) # formating data
#'     Loricariidae # informing Family of Curculionichthys insperatus
#'     Siluriformes # informing Order of Curculionichthys insperatus
#'     res_test <- whichFishAdd(data = taxon_data)
#' }
#' 
#' 
whichFishAdd <- function(data){
  tree_complete <- fishtree::fishtree_phylogeny()
  round_1_check <- match(data$s, tree_complete$tip.label)
  round_1_check <- round_1_check[!is.na(round_1_check)]
  if(length(round_1_check) == length(data$s)){
    nsertions <- rep("NA", nrow(data))
    data_insertions <- cbind(data, insertions)
    data_insertions[match(spp_on_tree, data$s), "insertions"] <- rep("Present_in_Tree", 
                                                                     nrow(data))
    return(data_insertions)
  }
  if(length(round_1_check) != length(data$s)){
    rank_order <- as.character(unique(data$o))
    rank_family <- as.character(unique(data$f))
    spp <- as.character(data$s)
    cichliformes_ord <- which(rank_order == "Cichliformes")
    if (length(cichliformes_ord) == 1) {
      rank_order[cichliformes_ord] <- "Perciformes"
    }
    all_families <- unique(unlist(lapply(rank_order, function(x) {
      fishbase[which(x == fishbase$Order), 10]
    })))
    families_in_orders <- suppressWarnings(all_families[which(unique(data$f) != 
                                                                all_families)])
    families_order_and_data <- unique(c(rank_family, families_in_orders))
    list_family <- vector(mode = "list", length = length(families_order_and_data))
    for (i in 1:length(families_order_and_data)) {
      list_family[[i]] <- tryCatch(paste(fishtree::fishtree_phylogeny(rank = families_order_and_data[i], 
                                                                            type = "chronogram_mrca")$tip.label), error = function(e) paste(print(families_order_and_data[i])))
    }
    names(list_family) <- families_order_and_data
    monotipic_family <- names(unlist(lapply(list_family, 
                                            function(x) which(length(x) == 1))))
    list_monotipic <- vector(mode = "list", length = length(monotipic_family))
    for (i in 1:length(monotipic_family)) {
      list_monotipic[[i]] <- tryCatch(fishtree::fishtree_taxonomy(rank = monotipic_family[i])[[1]]$taxonomy[[9]], 
                                      error = function(e) paste("not.found", "_", monotipic_family[i], 
                                                                sep = ""))
    }
    orders_to_add <- unique(unlist(list_monotipic[-which(sub("_.*", 
                                                             "", unlist(list_monotipic)) == "not.found")]))
    differences_orders_toadd <- setdiff(rank_order, orders_to_add)
    if (length(differences_orders_toadd) >= 1) {
      all_orders_include <- c(differences_orders_toadd, 
                              orders_to_add)
    }
    all_orders_include <- unique(c(rank_order, unique(orders_to_add)))
    list_order <- vector(mode = "list", length = length(all_orders_include))
    for (i in 1:length(all_orders_include)) {
      list_order[[i]] <- tryCatch(paste(fishtree::fishtree_phylogeny(rank = all_orders_include[i], 
                                                                           type = "chronogram_mrca")$tip.label), error = function(e) paste(print(all_orders_include[i])))
    }
    names(list_order) <- all_orders_include
    phylo_order <- filter_rank(order = list_order)
    phylo_order <- ape::makeNodeLabel(phy = phylo_order)
    order_rm_list <- names(unlist(lapply(list_order, function(x) which(length(x) == 
                                                                         1))))
    list_order <- list_order[-match(order_rm_list, names(list_order))]
    
    list_non_monotipic <- list_family[setdiff(names(list_family), 
                                              monotipic_family)]
    for (i in 1:length(list_non_monotipic)) {
      phylo_order <- ape::makeNodeLabel(phylo_order, "u", 
                                        nodeList = list(Fam_name = list_non_monotipic[[i]]))
      phylo_order$node.label[which(phylo_order$node.label == 
                                     "Fam_name")] <- paste(names(list_non_monotipic)[i])
    }
    families_in_tree <- families_order_and_data[which(!is.na(match(families_order_and_data, 
                                                                   phylo_order$node.label)) == T)]
    families_monotipic_notfound <- setdiff(monotipic_family, 
                                           families_in_tree)
    for (i in 1:length(families_monotipic_notfound)) {
      spp_tmp <- tryCatch(fishtree::fishtree_taxonomy(rank = families_monotipic_notfound[i])[[1]]$species, 
                          error = function(e) paste("not.found", "_", families_monotipic_notfound[i], 
                                                    sep = ""))
      spp_tmp <- gsub("\\ ", "_", spp_tmp)
      list_family[which(families_monotipic_notfound[i] == 
                          names(list_family))] <- list(spp_tmp)
    }
    phylo_order <- suppressWarnings(filter_rank(order = list_family))
    phylo_order <- ape::makeNodeLabel(phy = phylo_order)
    for (i in 1:length(list_order)) {
      phylo_order <- ape::makeNodeLabel(phylo_order, "u", 
                                        nodeList = list(Ord_name = list_order[[i]]))
      phylo_order$node.label[which(phylo_order$node.label == 
                                     "Ord_name") ] <- names(list_order)[i]
    }
    phylo_order <- phytools::force.ultrametric(phylo_order)
    families_not_found_fishtree <- names(unlist(lapply(lapply(list_family, 
                                                              function(x) {
                                                                sub("_.*", "", x)
                                                              }), function(y) which(length(y) == 1) & which(y == 
                                                                                                              "not.found"))))
    list_family_tobeaddnames <- list_family[-match(families_not_found_fishtree, 
                                                   names(list_family))]
    family_no_spp_in_tree <- names(unlist(lapply(lapply(list_family_tobeaddnames, 
                                                        function(x) {
                                                          sum(!is.na(match(x, phylo_order$tip.label)))
                                                        }), function(y) which(y == 0))))
    list_family_tobeaddnames <- list_family_tobeaddnames[-match(family_no_spp_in_tree, 
                                                                names(list_family_tobeaddnames))]
    for (i in 1:length(list_family_tobeaddnames)) {
      na_check <- sum(!is.na(match(list_family_tobeaddnames[[i]], 
                                   phylo_order$tip.label)))
      if (na_check == 1) {
        spp_singleton <- unlist(list(list_family_tobeaddnames[[i]][!is.na(match(list_family_tobeaddnames[[i]], 
                                                                                phylo_order$tip.label))]))
        spp_singleton_add <- paste(sub("_.*", "", spp_singleton), 
                                   "_", "singleton", sep = "")
        phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                      species = spp_singleton_add)
        list_family_tobeaddnames[i] <- list(c(spp_singleton, 
                                              spp_singleton_add))
      }
      phylo_order <- ape::makeNodeLabel(phylo_order, "u", 
                                        nodeList = list(Fam_name = list_family_tobeaddnames[[i]]))
      phylo_order$node.label[which(phylo_order$node.label == 
                                     "Fam_name") ] <- paste(names(list_family_tobeaddnames)[i])
    }
    spp_data <- 1:length(spp)
    names(spp_data) <- spp
    insert_spp <- treedata_modif(phy = phylo_order, data = spp_data, 
                                 warnings = F)$nc$data_not_tree
    if (length(insert_spp) >= 1) {
      genus_in_tree <- sub("_.*", "", phylo_order$tip.label)[match(sub("_.*", 
                                                                       "", insert_spp), sub("_.*", "", phylo_order$tip.label))][!is.na(sub("_.*", 
                                                                                                                                           "", phylo_order$tip.label)[match(sub("_.*", "", 
                                                                                                                                                                                insert_spp), sub("_.*", "", phylo_order$tip.label))])]
      species_to_genus1 <- insert_spp[which(is.na(insert_spp[match(sub("_.*", 
                                                                       "", insert_spp), genus_in_tree)]) == FALSE)]
      
      insert_spp2 <- setdiff(names(spp_data), c(phylo_order$tip.label, species_to_genus1))
      
      if (length(insert_spp2) == 0) {
        insertions <- rep("NA", nrow(data))
        data_insertions <- cbind(data, insertions)
        data_insertions[which(species_to_genus1 == data$s), 
                        "insertions"] <- rep("Congeneric_insertion", 
                                             length(species_to_genus1))
        spp_on_tree <- data[-match(species_to_genus1, 
                                   data$s), "s"]
        data_insertions[match(spp_on_tree, data$s), 
                        "insertions"] <- rep("Present_in_Tree", length(spp_on_tree))
        
        return(data_insertions)
      }
      if (length(insert_spp2) >= 1) {
        data_exRound2 <- data[match(insert_spp2, as.character(data$s)), 
        ]
        rank_family2 <- unique(as.character(data[match(insert_spp2, 
                                                       as.character(data$s)), 2]))
        list_spp_step2 <- vector(mode = "list", length = length(rank_family2))
        for (i in 1:length(rank_family2)) {
          list_spp_step2[[i]] <- tryCatch(paste(ape::extract.clade(phy = phylo_order, 
                                                                   node = as.character(rank_family2[i]))$tip.label), 
                                          error = function(e) paste("noFamily", as.character(data[which(rank_family2[i] == 
                                                                                                          data$f), 1]), sep = "_"))
        }
        names(list_spp_step2) <- rank_family2
        data_exRound3 <- data_exRound2[!is.na(match(data_exRound2$f, names(which(unlist(lapply(lapply(list_spp_step2, 
                                                                                                      function(x) which(sub("_.*", "", x) == "noFamily")), 
                                                                                               function(y) length(y))) > 0)))), ]
        spp_family <- 1:nrow(data_exRound2)
        names(spp_family) <- data_exRound2$s
        spp_with_family <- names(which(unlist(lapply(lapply(list_spp_step2, 
                                                            function(x) which(sub("_.*", "", x) != "noFamily")), 
                                                     function(y) which(length(y) != 0))) > 0))
        spp_family_inTree <- list_spp_step2[match(spp_with_family, names(list_spp_step2)
        )]
        spp_to_add_round2 <- setdiff(data_exRound2$s, 
                                     data_exRound3$s)
        
        if (dim(data_exRound3)[1] == 0) {
          insertions <- rep("NA", nrow(data))
          data_insertions <- cbind(data, insertions)
          data_insertions[match(species_to_genus1, data$s), "insertions"] <- rep("Congeneric_insertion", length(species_to_genus1))
          data_insertions[match(data_exRound2$s, data$s), 
                          "insertions"] <- rep("Family_insertion", 
                                               length(data_exRound2$s))
          spp_on_tree <- data[-match(c(species_to_genus1,
                                       data$s[match(data_exRound2$s, data$s)]), 
                                     data$s), "s"]
          data_insertions[match(spp_on_tree, data$s), 
                          "insertions"] <- rep("Present_in_Tree", 
                                               length(spp_on_tree))
          
          return(data_insertions)
        }
        else {
          insertions <- rep("NA", nrow(data))
          data_insertions <- cbind(data, insertions)
          species_to_genus2 <- unlist(species_to_genus2)
          data_insertions <- cbind(data, insertions)
          data_insertions[match(species_to_genus1, data$s), "insertions"] <- rep("Congeneric_insertion", length(species_to_genus1))
          if(length(species_to_genus2) >= 1){
            data_insertions[match(species_to_genus2, data$s), "insertions"] <- rep("Congeneric_insertion_roundFamily", length(unlist(species_to_genus2)))
            data_exRound2 <- data_exRound2[-match(species_to_genus2, data_exRound2$s), ]
          }
          family_level_insertions <- unique(setdiff(data_exRound2$s, 
                                                    data_exRound3$s))
          family_insertions <- setdiff(family_level_insertions, species_to_genus1)
          data_insertions[match(family_insertions, 
                                data$s), "insertions"] <- rep("Family_insertion", 
                                                              length(family_insertions))            
          
          data_insertions[match(data_exRound3$s, data$s), 
                          "insertions"] <- rep("Order_insertion", 
                                               length(data_exRound3$s))
          spp_on_tree <- data[-match(c(species_to_genus1, species_to_genus2,
                                       family_insertions, 
                                       data_exRound3$s), data$s), "s"]
          data_insertions[match(spp_on_tree, data$s), 
                          "insertions"] <- rep("Present_in_Tree", 
                                               length(spp_on_tree))
          
          return(data_insertions)
        }
      }
    }
  }
}

