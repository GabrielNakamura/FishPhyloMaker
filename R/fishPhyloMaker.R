#' Obtaining fish phylogeny according to a local pool of species
#'
#' @param data A data frame with three columns containing the name of species (s), the Family (f) and the Order (o). This data frame can be generated 
#'     with tab_function function.
#' @param return.insertions Logical, if TRUE (default) the output is a list of length two containing the phylogeny 
#'     and a dataframe with a column indicating at which level each species was inserted.
#' @param insert.base.node Logical argument indicating if the species must be added automatically
#'     in the family and order (when needed) nodes. Default is FALSE
#' @param progress.bar Logical argument. If TRUE (default) a progress bar will be shown in console.
#'    
#'     
#' @return A newick object containing the phylogeny with the species in data object. If return.insertions = TRUE the output
#'     will be a list of length two containing the newick phylogeny and a data frame equal that provided in data plus a column 
#'     indicating at which level each species was inserted in the tree.
#' @export
#'
#' @examples \donttest{
#'     data("taxon_data_PhyloMaker")
#'     res_phylo <- FishPhyloMaker(data = taxon_data_PhyloMaker,
#'     insert.base.node = TRUE, 
#'     return.insertions = TRUE, 
#'     progress.bar = TRUE)
#' }
#'     
#' 

FishPhyloMaker <- function (data, 
                            insert.base.node = FALSE, 
                            return.insertions = TRUE, 
                            progress.bar = TRUE) 
{
  if (dim(data)[2] != 3) {
    stop("/n data must be a dataframe with three columns (s, f, o)")
  }
  if (any(is.na(match(colnames(data), c("s", "f", "o"))))) {
    stop("/n Columns of data object must be named with s, f, and o letters")
  }
  if (is.data.frame(data) == FALSE) {
    stop("/n data must be a data frame object")
  }
  fishbasedata <- as.data.frame(data.frame(rfishbase::load_taxa()))
  tree_complete <- fishtree::fishtree_phylogeny()
  round_1_check <- match(data$s, tree_complete$tip.label)
  round_1_check <- round_1_check[!is.na(round_1_check)]
  if (length(round_1_check) == length(data$s)) {
    data_final <- 1:length(as.character(data$s))
    names(data_final) <- as.character(data$s)
    tree_res <- suppressWarnings(ape::drop.tip(phy = tree_complete, 
                                               tip = treedata_modif(phy = tree_complete, data = data_final)$nc$tree_not_data))
    if (return.insertions == TRUE) {
      insertions <- rep("NA", nrow(data))
      data_insertions <- cbind(data, insertions)
      data_insertions[, "insertions"] <- rep("Present_in_Tree", 
                                             nrow(data))
      list_res <- vector(mode = "list", length = 2)
      list_res[[1]] <- tree_res
      list_res[[2]] <- data_insertions
      names(list_res) <- c("Phylogeny", "Insertions_data")
      return(list_res)
    }
    else {
      return(tree_res)
    }
  }
  if (length(round_1_check) != length(data$s)) {
    rank_order <- as.character(unique(data$o))
    rank_family <- as.character(unique(data$f))
    spp <- as.character(data$s)
    all_families <- unique(unlist(lapply(rank_order, function(x) {
      fishbasedata[which(x == fishbasedata$Order), 5]
    })))
    families_in_orders <- suppressWarnings(all_families[which(unique(data$f) != 
                                                                all_families)])
    families_order_and_data <- unique(c(rank_family, families_in_orders))
    list_family <- vector(mode = "list", length = length(families_order_and_data))
    if (progress.bar == TRUE) {
      pb_download_family <- progress::progress_bar$new(format = "  downloading families [:bar] :percent", 
                                                       total = length(families_order_and_data), clear = FALSE, 
                                                       width = 60, current = "<", incomplete = ">", 
                                                       complete = ">")
    }
    for (i in 1:length(families_order_and_data)) {
      list_family[[i]] <- tryCatch(paste(fishtree::fishtree_phylogeny(rank = families_order_and_data[i], 
                                                                      type = "chronogram_mrca")$tip.label), error = function(e) paste(families_order_and_data[i]))
      if (progress.bar == TRUE) {
        pb_download_family$tick()
      }
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
    if (progress.bar == TRUE) {
      pb_download_order <- progress::progress_bar$new(format = "  downloading orders [:bar] :percent", 
                                                      total = length(all_orders_include), clear = FALSE, 
                                                      width = 60, current = "<", incomplete = ">", 
                                                      complete = ">")
    }
    for (i in 1:length(all_orders_include)) {
      list_order[[i]] <- tryCatch(paste(fishtree::fishtree_phylogeny(rank = all_orders_include[i], 
                                                                     type = "chronogram_mrca")$tip.label), error = function(e) paste(print(all_orders_include[i])))
      if (progress.bar == TRUE) {
        pb_download_order$tick()
      }
    }
    names(list_order) <- all_orders_include
    phylo_order <- filter_rank(order = list_order)
    phylo_order <- ape::makeNodeLabel(phy = phylo_order)
    order_rm_list <- names(unlist(lapply(list_order, function(x) which(length(x) == 
                                                                         1))))
    list_order <- list_order[-match(order_rm_list, names(list_order))]
    list_non_monotipic <- list_family[setdiff(names(list_family), 
                                              monotipic_family)]
    if (progress.bar == TRUE) {
      pb_names_family <- progress::progress_bar$new(format = "  Naming family nodes [:bar] :percent", 
                                                    total = length(list_non_monotipic), clear = FALSE, 
                                                    width = 60, current = "<", incomplete = ">", 
                                                    complete = ">")
    }
    for (i in 1:length(list_non_monotipic)) {
      phylo_order <- ape::makeNodeLabel(phylo_order, "u", 
                                        nodeList = list(Fam_name = list_non_monotipic[[i]]))
      phylo_order$node.label[which(phylo_order$node.label == 
                                     "Fam_name")] <- paste(names(list_non_monotipic)[i])
      if (progress.bar == TRUE) {
        pb_names_family$tick()
      }
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
    if (length(list_family_tobeaddnames) > 0) {
      if (progress.bar == TRUE) {
        pb_names_family_toadd <- progress::progress_bar$new(format = "  Naming monotipic family nodes [:bar] :percent", 
                                                            total = length(list_family_tobeaddnames), clear = FALSE, 
                                                            width = 60, current = "<", incomplete = ">", 
                                                            complete = ">")
      }
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
        phylo_order <- ape::makeNodeLabel(phylo_order, 
                                          "u", nodeList = list(Fam_name = list_family_tobeaddnames[[i]]))
        phylo_order$node.label[which(phylo_order$node.label == 
                                       "Fam_name")] <- paste(names(list_family_tobeaddnames)[i])
        if (progress.bar == TRUE) {
          pb_names_family_toadd$tick()
        }
      }
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
      if (length(species_to_genus1) >= 1) {
        if (progress.bar == TRUE) {
          pb_congeneric <- progress::progress_bar$new(format = "Adding congeneric species [:bar] :percent", 
                                                      total = length(species_to_genus1), clear = FALSE, 
                                                      width = 60, current = "<", incomplete = ">", 
                                                      complete = ">")
        }
        for (i in 1:length(species_to_genus1)) {
          phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                        species = species_to_genus1[i])
          if (progress.bar == TRUE) {
            pb_congeneric$tick()
          }
        }
      }
      insert_spp2 <- treedata_modif(phy = phylo_order, 
                                    data = spp_data, warnings = F)$nc$data_not_tree
      if (length(insert_spp2) == 0) {
        data_final <- 1:length(as.character(data$s))
        names(data_final) <- as.character(data$s)
        tree_res <- suppressWarnings(ape::drop.tip(phy = phylo_order, 
                                                   tip = treedata_modif(phy = phylo_order, data = data_final)$nc$tree_not_data))
        if (return.insertions == TRUE) {
          insertions <- rep("NA", nrow(data))
          data_insertions <- cbind(data, insertions)
          data_insertions[match(species_to_genus1, data$s), 
                          "insertions"] <- rep("Congeneric_insertion", 
                                               length(species_to_genus1))
          spp_on_tree <- data[-match(species_to_genus1, 
                                     data$s), "s"]
          data_insertions[match(spp_on_tree, data$s), 
                          "insertions"] <- rep("Present_in_Tree", length(spp_on_tree))
          list_res <- vector(mode = "list", length = 2)
          list_res[[1]] <- tree_res
          list_res[[2]] <- data_insertions
          names(list_res) <- c("Phylogeny", "Insertions_data")
          return(list_res)
        }
        else {
          return(tree_res)
        }
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
        data_exRound3 <- data_exRound2[!is.na(match(data_exRound2$f, 
                                                    names(which(unlist(lapply(lapply(list_spp_step2, 
                                                                                     function(x) which(sub("_.*", "", x) == "noFamily")), 
                                                                              function(y) length(y))) > 0)))), ]
        spp_family <- 1:nrow(data_exRound2)
        names(spp_family) <- data_exRound2$s
        spp_with_family <- names(which(unlist(lapply(lapply(list_spp_step2, 
                                                            function(x) which(sub("_.*", "", x) != "noFamily")), 
                                                     function(y) which(length(y) != 0))) > 0))
        spp_family_inTree <- list_spp_step2[match(spp_with_family, 
                                                  names(list_spp_step2))]
        spp_to_add_round2 <- setdiff(data_exRound2$s, 
                                     data_exRound3$s)
        if (insert.base.node == TRUE) {
          pb_length <- unique(sub("_.*", "", as.character(spp_to_add_round2)))
          if (progress.bar == TRUE) {
            pb_insert_family_node <- progress::progress_bar$new(format = "Adding species to family node [:bar] :percent", 
                                                                total = length(pb_length), clear = FALSE, 
                                                                width = 60, current = "<", incomplete = ">", 
                                                                complete = ">")
          }
          count <- 0
          species_to_genus2 <- vector(mode = "list")
          while (length(spp_to_add_round2) >= 1) {
            count <- count + 1
            family_name <- data[match(spp_to_add_round2[1], 
                                      data$s), "f"]
            node_family <- which(c(phylo_order$tip.label, 
                                   phylo_order$node.label) == family_name)
            phylo_order <- phytools::bind.tip(tree = phylo_order, 
                                              tip.label = spp_to_add_round2[1], where = node_family, 
                                              position = 0)
            spp_data <- 1:length(spp_to_add_round2)
            names(spp_data) <- spp_to_add_round2
            insert_spp <- treedata_modif(phy = phylo_order, 
                                         data = spp_data, warnings = F)$nc$data_not_tree
            genus_in_tree <- sub("_.*", "", phylo_order$tip.label)[match(sub("_.*", 
                                                                             "", insert_spp), sub("_.*", "", phylo_order$tip.label))][!is.na(sub("_.*", 
                                                                                                                                                 "", phylo_order$tip.label)[match(sub("_.*", 
                                                                                                                                                                                      "", insert_spp), sub("_.*", "", phylo_order$tip.label))])]
            if (length(genus_in_tree) >= 1) {
              species_to_genus <- insert_spp[which(is.na(insert_spp[match(sub("_.*", 
                                                                              "", insert_spp), genus_in_tree)]) == 
                                                     FALSE)]
              species_to_genus2[[count]] <- species_to_genus
              for (i in 1:length(species_to_genus)) {
                phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                              species = species_to_genus[i])
              }
              insert_spp <- treedata_modif(phy = phylo_order, 
                                           data = spp_data, warnings = F)$nc$data_not_tree
            }
            spp_to_add_round2 <- setdiff(insert_spp, 
                                         data_exRound3$s)
            if (progress.bar == TRUE) {
              pb_insert_family_node$tick()
            }
          }
        }
        else {
          count <- 0
          species_to_genus2 <- vector(mode = "list")
          while (length(spp_to_add_round2) >= 1) {
            count <- count + 1
            family_name <- data[match(spp_to_add_round2[1], 
                                      data$s), "f"]
            spp_Byfamily_inTree <- as.character(unlist(spp_family_inTree[match(family_name, 
                                                                               names(spp_family_inTree))]))
            user_option_spp <- unique(sub("_.*", "", 
                                          as.character(unlist(spp_Byfamily_inTree))))
            local_to_add_spp <- readline(prompt = print_cat(print_cat = user_option_spp, 
                                                            spp = spp_to_add_round2[1], family_name))
            spp_user_opt <- unlist(strsplit(local_to_add_spp, 
                                            split = " "))
            if (length(spp_user_opt) > 2) {
              stop("/n Only two genus can be chosen to insert species")
            }
            if (length(spp_user_opt) < 1) {
              stop("/n At least one genus/family can be chosen to insert species")
            }
            if (any(is.na(match(spp_user_opt, user_option_spp) && 
                          is.na(match(spp_user_opt, family_name))))) {
              stop("/n Choose a validy genus/family to insert species")
            }
            if (length(spp_user_opt) == 2) {
              spp_to_add_tmp <- spp_Byfamily_inTree[match(spp_user_opt, 
                                                          sub("_.*", "", as.character(unlist(spp_Byfamily_inTree))))]
              node_btw_genus <- phytools::fastMRCA(tree = phylo_order, 
                                                   sp1 = spp_to_add_tmp[1], sp2 = spp_to_add_tmp[2])
              phylo_order <- phytools::bind.tip(tree = phylo_order, 
                                                tip.label = spp_to_add_round2[1], where = node_btw_genus, 
                                                position = 0)
            }
            if (length(spp_user_opt) == 1) {
              if (length(match(spp_user_opt, c(user_option_spp, 
                                               family_name))) != 1) {
                stop(paste("\n Check the spelling of Genus or Family in:", 
                           spp_user_opt, sep = " "))
              }
              if (any(spp_user_opt == family_name)) {
                node_family <- which(c(phylo_order$tip.label, 
                                       phylo_order$node.label) == family_name)
                phylo_order <- phytools::bind.tip(tree = phylo_order, 
                                                  tip.label = spp_to_add_round2[1], where = node_family, 
                                                  position = 0)
              }
              else {
                genus_nspp <- match(spp_user_opt, names(which(table(sub("_.*", 
                                                                        "", as.character(unlist(spp_Byfamily_inTree)))) > 
                                                                1)))
                if (!is.na(genus_nspp) == TRUE) {
                  list_to_be_named <- list(spp_user_opt = spp_user_opt)
                  name_list <- paste("MRCA", count, sep = "")
                  names(list_to_be_named) <- name_list
                  phylo_order2 <- ape::makeNodeLabel(phy = phylo_order, 
                                                     method = "u", nodeList = list_to_be_named)
                  position_MRCA2 <- which(c(phylo_order2$tip.label, 
                                            phylo_order2$node.label) == name_list)
                  check_name_node <- phylo_order$node.label[position_MRCA2 - 
                                                              length(phylo_order$tip.label)]
                  if (check_name_node == family_name) {
                    list_to_be_named <- list(spp_user_opt = spp_user_opt)
                    name_list <- paste("MRCA", count, 
                                       sep = "")
                    names(list_to_be_named) <- name_list
                    phylo_order <- ape::makeNodeLabel(phy = phylo_order, 
                                                      method = "u", nodeList = list_to_be_named)
                    position_MRCA2 <- which(c(phylo_order$tip.label, 
                                              phylo_order$node.label) == name_list)
                    size_branch <- phylo_order$edge.length[sapply(position_MRCA2, 
                                                                  function(x, y) which(y == x), y = phylo_order$edge[, 
                                                                                                                     2])]
                    phylo_order <- phytools::bind.tip(phylo_order, 
                                                      spp_to_add_round2[1], where = position_MRCA2, 
                                                      position = size_branch/2)
                    node_pos_new <- phytools::fastMRCA(tree = phylo_order, 
                                                       sp1 = spp_Byfamily_inTree[1], sp2 = spp_to_add_round2[1])
                    phylo_order$node.label[node_pos_new - 
                                             length(phylo_order$tip.label)] <- family_name
                  }
                  else {
                    phylo_order <- ape::makeNodeLabel(phy = phylo_order, 
                                                      method = "u", nodeList = list_to_be_named)
                    position_MRCA <- which(c(phylo_order$tip.label, 
                                             phylo_order$node.label) == name_list)
                    size_branch <- phylo_order$edge.length[sapply(position_MRCA, 
                                                                  function(x, y) which(y == x), y = phylo_order$edge[, 
                                                                                                                     2])]
                    phylo_order <- phytools::bind.tip(phylo_order, 
                                                      spp_to_add_round2[1], where = position_MRCA, 
                                                      position = size_branch/2)
                  }
                }
                if (is.na(genus_nspp) == TRUE) {
                  phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                                species = paste(sub("_.*", "", as.character(spp_user_opt)), 
                                                                                "toadd", sep = "_"))
                  position_problem2 <- which(phylo_order$tip.label == 
                                               paste(sub("_.*", "", as.character(spp_user_opt)), 
                                                     "toadd", sep = "_"))
                  phylo_order$tip.label[position_problem2] <- spp_to_add_round2[1]
                }
              }
            }
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
                                                                              "", insert_spp), genus_in_tree)]) == 
                                                     FALSE)]
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
              break
            }
            else {
              list_spp_step2 <- vector(mode = "list", 
                                       length = length(rank_family2))
              for (i in 1:length(rank_family2)) {
                list_spp_step2[[i]] <- tryCatch(paste(ape::extract.clade(phy = phylo_order, 
                                                                         node = as.character(rank_family2[i]))$tip.label), 
                                                error = function(e) paste("noFamily", 
                                                                          as.character(data[which(rank_family2[i] == 
                                                                                                    data$f), 1]), sep = "_"))
              }
              names(list_spp_step2) <- rank_family2
              spp_family <- 1:nrow(data_exRound2)
              names(spp_family) <- data_exRound2$s
              spp_with_family <- names(which(unlist(lapply(lapply(list_spp_step2, 
                                                                  function(x) which(sub("_.*", "", x) != 
                                                                                      "noFamily")), function(y) which(length(y) != 
                                                                                                                        0))) > 0))
              spp_family_inTree <- list_spp_step2[match(spp_with_family, 
                                                        names(list_spp_step2))]
              spp_to_add_round2 <- setdiff(insert_spp, 
                                           data_exRound3$s)
            }
          }
        }
        if (is.null(data_exRound3) | dim(data_exRound3)[1] == 0) {
          if (return.insertions == TRUE) {
            data_final <- 1:length(as.character(data$s))
            names(data_final) <- as.character(data$s)
            tree_res <- suppressWarnings(ape::drop.tip(phy = phylo_order, 
                                                       tip = treedata_modif(phy = phylo_order, 
                                                                            data = data_final)$nc$tree_not_data))
            insertions <- rep("NA", nrow(data))
            data_insertions <- cbind(data, insertions)
            Present_in_tree <- data$s[!is.na(match(data$s, 
                                                   tree_complete$tip.label))]
            Congeneric_insertion <- species_to_genus1
            not_inserted <- data$s[is.na(match(data$s, 
                                               tree_res$tip.label))]
            Congeneric_round_family <- unlist(species_to_genus2)
            if(is.null(Congeneric_round_family) & length(not_inserted) == 0){
              data_insertions[match(insert_spp2, data_insertions$s), "insertions"] <- "Family_insertion"
            } else {
              Family_insertion <- insert_spp2[-match(c(Congeneric_round_family, 
                                                       not_inserted), insert_spp2)]
              data_insertions[match(Family_insertion, data_insertions$s), 
                              "insertions"] <- "Family_insertion"
            }
            
            data_insertions[match(Present_in_tree, data_insertions$s), 
                            "insertions"] <- "Present_in_Tree"
            data_insertions[match(Congeneric_insertion, 
                                  data_insertions$s), "insertions"] <- "Congeneric_insertion"
            
            if (length(not_inserted) >= 1) {
              data_insertions[match(not_inserted, data_insertions$s), 
                              "insertions"] <- "Not_inserted"
            }
            if (length(Congeneric_round_family) >= 1) {
              data_insertions[match(Congeneric_round_family, 
                                    data_insertions$s), "insertions"] <- "Congeneric_Family_level"
            }
            list_res <- vector(mode = "list", length = 2)
            list_res[[1]] <- tree_res
            list_res[[2]] <- data_insertions
            names(list_res) <- c("Phylogeny", "Insertions_data")
            return(list_res)
          }
          else {
            return(tree_res)
          }
        }
        else {
          data_exRound3 <- data_exRound3[is.na(match(data_exRound3$o, 
                                                     order_rm_list)), ]
          rank_order_Round3 <- rank_order[match(data_exRound3$o, 
                                                rank_order)]
          families_round3 <- lapply(lapply(rank_order_Round3, 
                                           function(x) {
                                             fishbasedata[which(x == fishbasedata$Order), 
                                                          5]
                                           }), function(y) unique(y))
          names(families_round3) <- rank_order_Round3
          orders_round_3 <- unique(rank_order_Round3)
          families_round3_check <- unique(as.character(unlist(families_round3)))
          check_round3_insertion <- match(families_round3_check, 
                                          families_in_tree)
          check_allNA_round3 <- unique(check_round3_insertion)
          if (all(is.na(check_allNA_round3)) == TRUE) {
            no_represent_tree <- data_exRound3$s
            data_final <- 1:length(as.character(data$s))
            names(data_final) <- as.character(data$s)
            tree_res <- suppressWarnings(ape::drop.tip(phy = phylo_order, 
                                                       tip = treedata_modif(phy = phylo_order, 
                                                                            data = data_final)$nc$tree_not_data))
            if (return.insertions == TRUE) {
              insertions <- rep("NA", nrow(data))
              species_to_genus2 <- unlist(species_to_genus2)
              data_insertions <- cbind(data, insertions)
              Present_in_tree <- data$s[!is.na(match(data$s, 
                                                     tree_complete$tip.label))]
              Congeneric_insertion <- species_to_genus1
              not_inserted <- data$s[is.na(match(data$s, 
                                                 tree_res$tip.label))]
              Congeneric_round_family <- species_to_genus2
              if(is.null(Congeneric_round_family) & length(not_inserted) == 0){
                data_insertions[match(insert_spp2, data_insertions$s), "insertions"] <- "Family_insertion"
              } else {
                Family_insertion <- insert_spp2[-match(c(Congeneric_round_family, 
                                                         not_inserted), insert_spp2)]
                data_insertions[match(Family_insertion, data_insertions$s), 
                                "insertions"] <- "Family_insertion"
              }
              data_insertions[match(Present_in_tree, 
                                    data_insertions$s), "insertions"] <- "Present_in_Tree"
              data_insertions[match(Congeneric_insertion, 
                                    data_insertions$s), "insertions"] <- "Congeneric_insertion"
              if (length(not_inserted) >= 1) {
                data_insertions[match(not_inserted, data_insertions$s), 
                                "insertions"] <- "Not_inserted"
              }
              if (length(Congeneric_round_family) >= 
                  1) {
                data_insertions[match(Congeneric_round_family, 
                                      data_insertions$s), "insertions"] <- "Congeneric_Family_level"
              }
              list_res <- vector(mode = "list", length = 2)
              list_res[[1]] <- tree_res
              list_res[[2]] <- data_insertions
              names(list_res) <- c("Phylogeny", "Insertions_data")
              return(list_res)
            }
            return(tree_res)
          }
          order_add_round3 <- unique(data_exRound3$o)
          if (progress.bar == TRUE) {
            pb_insert_order_node_remaining <- progress::progress_bar$new(format = "Adding species to order node [:bar] :percent", 
                                                                         total = nrow(data_exRound3), clear = FALSE, 
                                                                         width = 60, current = "<", incomplete = ">", 
                                                                         complete = ">")
          }
          for (i in 1:length(order_add_round3)) {
            list_names_round3 <- list_order[[which(names(list_order) == 
                                                     order_add_round3[i])]]
            phylo_order <- ape::makeNodeLabel(phylo_order, 
                                              "u", nodeList = list(Ord_name = list_names_round3))
            phylo_order$node.label[which(phylo_order$node.label == 
                                           "Ord_name")] <- order_add_round3[i]
          }
          if (insert.base.node == TRUE) {
            for (i in 1:nrow(data_exRound3)) {
              node_order_pos <- which(c(phylo_order$tip.label, 
                                        phylo_order$node.label) == data_exRound3$o[i])
              phylo_order <- phytools::bind.tip(tree = phylo_order, 
                                                tip.label = data_exRound3$s[i], where = node_order_pos, 
                                                position = 0)
              pb_insert_order_node_remaining$tick()
            }
          }
          else {
            for (i in 1:length(families_round3)) {
              user_option_family <- phylo_order$node.label[match(unlist(families_round3[[i]]), 
                                                                 phylo_order$node.label)[-which(match(match(unlist(families_round3[[i]]), 
                                                                                                            phylo_order$node.label), NA) == 1)]]
              local_to_add_spp_family <- readline(prompt = print_cat_family(print_cat = unlist(user_option_family), 
                                                                            spp = data_exRound3$s[i], data_exRound3$o[i]))
              family_user_opt <- unlist(strsplit(local_to_add_spp_family, 
                                                 split = " "))
              if (length(family_user_opt) == 1) {
                if (family_user_opt == data_exRound3$o[i]) {
                  node_order_pos <- which(c(phylo_order$tip.label, 
                                            phylo_order$node.label) == data_exRound3$o[i])
                  phylo_order <- phytools::bind.tip(tree = phylo_order, 
                                                    tip.label = data_exRound3$s[i], where = node_order_pos, 
                                                    position = 0)
                }
                family_nspp <- length(list_family[[match(family_user_opt, 
                                                         names(list_family))]])
                if (family_nspp > 1) {
                  position_family <- which(family_user_opt == 
                                             c(phylo_order$tip.label, phylo_order$node.label))
                  size_branch_family <- phylo_order$edge.length[sapply(position_family, 
                                                                       function(x, y) which(y == x), y = phylo_order$edge[, 
                                                                                                                          2])]
                  phylo_order <- phytools::bind.tip(phylo_order, 
                                                    data_exRound3$s[i], where = position_family, 
                                                    position = size_branch_family/2)
                }
                if (family_nspp == 1) {
                  phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                                species = paste(sub("_.*", "", as.character(list_family[[which(names(list_family) == 
                                                                                                                                 family_user_opt)]][1])), "toadd", 
                                                                                sep = "_"))
                  position_problem_family <- which(phylo_order$tip.label == 
                                                     paste(sub("_.*", "", as.character(list_family[[which(names(list_family) == 
                                                                                                            family_user_opt)]][1])), "toadd", 
                                                           sep = "_"))
                  phylo_order$tip.label[position_problem_family] <- data_exRound3$s[i]
                }
              }
              if (length(family_user_opt) == 2) {
                spp_to_add_tmp_family1 <- ape::extract.clade(phylo_order, 
                                                             node = family_user_opt[1])$tip.label[1]
                spp_to_add_tmp_family2 <- ape::extract.clade(phylo_order, 
                                                             node = family_user_opt[2])$tip.label[1]
                node_btw_genus_family <- phytools::fastMRCA(tree = phylo_order, 
                                                            sp1 = spp_to_add_tmp_family1, sp2 = spp_to_add_tmp_family2)
                phylo_order <- phytools::bind.tip(tree = phylo_order, 
                                                  tip.label = data_exRound3$s[i], where = node_btw_genus_family, 
                                                  position = 0)
              }
            }
          }
          data_final <- 1:length(as.character(data$s))
          names(data_final) <- as.character(data$s)
          tree_res <- suppressWarnings(ape::drop.tip(phy = phylo_order, 
                                                     tip = treedata_modif(phy = phylo_order, data = data_final)$nc$tree_not_data))
          if (return.insertions == TRUE) {
            insertions <- rep("NA", nrow(data))
            data_insertions <- cbind(data, insertions)
            species_to_genus2 <- unlist(species_to_genus2)
            data_insertions <- cbind(data, insertions)
            Present_in_tree <- data$s[!is.na(match(data$s, 
                                                   tree_complete$tip.label))]
            Congeneric_insertion <- species_to_genus1
            not_inserted <- data$s[is.na(match(data$s, 
                                               tree_res$tip.label))]
            Congeneric_round_family <- species_to_genus2
            if(is.null(Congeneric_round_family) & length(not_inserted) == 0){
              data_insertions[match(insert_spp2, data_insertions$s), "insertions"] <- "Family_insertion"
            } else {
              Family_insertion <- insert_spp2[-match(c(Congeneric_round_family, 
                                                       not_inserted), insert_spp2, nomatch = 0)]
              data_insertions[match(Family_insertion, data_insertions$s), 
                              "insertions"] <- "Family_insertion"
            }
            data_insertions[match(Present_in_tree, data_insertions$s), 
                            "insertions"] <- "Present_in_Tree"
            data_insertions[match(Congeneric_insertion, 
                                  data_insertions$s), "insertions"] <- "Congeneric_insertion"
            if (length(not_inserted) >= 1) {
              data_insertions[match(not_inserted, data_insertions$s), 
                              "insertions"] <- "Not_inserted"
            }
            if (length(Congeneric_round_family) >= 1) {
              data_insertions[match(Congeneric_round_family, 
                                    data_insertions$s), "insertions"] <- "Congeneric_Family_level"
            }
            if (dim(data_exRound3)[1] >= 1) {
              data_insertions[match(data_exRound3$s, 
                                    data_insertions$s), "insertions"] <- "Order_insertion"
            }
            list_res <- vector(mode = "list", length = 2)
            list_res[[1]] <- tree_res
            list_res[[2]] <- data_insertions
            names(list_res) <- c("Phylogeny", "Insertions_data")
            return(list_res)
          }
          else {
            return(tree_res)
          }
        }
      }
    }
  }
}
