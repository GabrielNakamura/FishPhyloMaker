
name_tree_nodes <- 
function(data){
  spp <- as.character(data$s)
  data_order_family <- get_phylo_order(data) # get species from all family and order  in the tree
  list_family <- data_order_family$family # sampled species in each family
  families_order_and_data <- data_order_family$families_order_and_data
  monotipic_family <- names(unlist(lapply(list_family, 
                                          function(x) which(length(x) == 1))))
  list_order <- data_order_family$order # sampled species in each order
  phylo_order <- filter_rank(order = list_order)
  phylo_order <- ape::makeNodeLabel(phy = phylo_order)
  order_rm_list <- names(unlist(lapply(list_order, function(x) which(length(x) == 
                                                                       1)))) # order with only one family to be removed
  list_order <- list_order[-match(order_rm_list, names(list_order))]
  list_non_monotipic <- list_family[setdiff(names(list_family), 
                                            monotipic_family)]
  
  for (i in 1:length(list_non_monotipic)) { # naming families with more than one species in tree
    node_anc <- ape::getMRCA(phy = phylo_order, tip = list_non_monotipic[[i]])
    phylo_order$node.label[node_anc - ape::Ntip(phylo_order)] <- names(list_non_monotipic)[i]
  }
  families_in_tree <- families_order_and_data[which(!is.na(match(families_order_and_data, 
                                                                 phylo_order$node.label)) == T)]
  families_monotipic_notfound <- setdiff(monotipic_family, 
                                         families_in_tree)
  for (i in 1:length(families_monotipic_notfound)) { # downloading species from monotipic families
    # i = 6
    spp_tmp <- tryCatch(fishtree::fishtree_taxonomy(rank = families_monotipic_notfound[i])[[1]]$species, 
                        error = function(e) paste("not.found", "_", families_monotipic_notfound[i], 
                                                  sep = ""))
    spp_tmp <- gsub("\\ ", "_", spp_tmp)
    list_family[which(families_monotipic_notfound[i] == 
                        names(list_family))] <- list(spp_tmp)
  }
  phylo_order <- suppressWarnings(filter_rank(order = list_family))
  phylo_order <- ape::makeNodeLabel(phy = phylo_order)
  phylo_order <- phytools::force.ultrametric(phylo_order)
  genus_monotip <- 
    lapply(list_family, function(x){
      sub("_.*", "", x)
    }
    )
  genus_not_found <- lapply(genus_monotip, function(y) which(length(y) == 1) & which(y == "not.found"))
  families_not_found_fishtree <- names(unlist(genus_not_found)) # names not in fishtree 
  list_family_tobeaddnames <- list_family[-match(families_not_found_fishtree, 
                                                 names(list_family))] # only families with at least one representative
  family_no_spp_tree <- lapply(list_family_tobeaddnames, 
                               function(x) {
                                 sum(!is.na(match(x, phylo_order$tip.label)))
                               })
  family_no_spp_in_tree <- names(unlist(lapply(family_no_spp_tree, function(y) which(y == 0))))
  list_family_tobeaddnames <- list_family_tobeaddnames[-match(family_no_spp_in_tree, 
                                                              names(list_family_tobeaddnames))]
  if (length(list_family_tobeaddnames) > 0) { # families to be added including monotipic families
    
    for (i in 1:length(list_family_tobeaddnames)) {
      na_check <- sum(!is.na(match(list_family_tobeaddnames[[i]], 
                                   phylo_order$tip.label)))
      if (na_check == 1) { # monotipic families
        spp_singleton <- unlist(list(list_family_tobeaddnames[[i]][!is.na(match(list_family_tobeaddnames[[i]], 
                                                                                phylo_order$tip.label))]))
        spp_singleton_add <- paste(sub("_.*", "", spp_singleton), 
                                   "_", "singleton", sep = "")
        phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                      species = spp_singleton_add)
        list_family_tobeaddnames[i] <- list(c(spp_singleton, 
                                              spp_singleton_add))
      }
      tip_spp <- list_family_tobeaddnames[[i]][!is.na(match(list_family_tobeaddnames[[i]], phylo_order$tip.label))]
      node_anc <- ape::getMRCA(phy = phylo_order, tip = tip_spp)
      phylo_order$node.label[node_anc - ape::Ntip(phylo_order)] <- names(list_family_tobeaddnames)[i]
    }
  }
}

