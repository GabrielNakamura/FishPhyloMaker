
#' Editing fish phylogeny according to a local pool of species
#'
#' @param data A data frame with three columns containing the name of species (s), the Family (f) and the Order (o). This data frame can be generated 
#'     with tab_function function
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
phyloMatch<- function(data){
  
  #organizing taxonomic levels
  rank_order <- as.character(unique(data$o)) #vector with orders
  rank_family <- as.character(unique(data$f)) #vector with families
  spp <- as.character(data$s) #vector with species
  list_family <- vector(mode = "list", length= length(rank_family))
  
  # list of families within genus presented in data
  cichliformes_ord <- which(rank_order == "Cichliformes")
  if(length(cichliformes_ord) == 1){
    rank_order[cichliformes_ord] <- "Perciformes"
  }
  all_families <- unique(unlist(lapply(rank_order, function(x){
    fishbase[which(x == fishbase$Order), 10]
  }))) # all families in data
  families_in_orders <- suppressWarnings(all_families[which(unique(data$f) != all_families)]) # all families which the orders are in data
  families_order_and_data <- c(rank_family, families_in_orders)
  
  #filtering for family
  for(i in 1:length(families_order_and_data)){
    list_family[[i]]<- tryCatch(paste(print(fishtree::fishtree_phylogeny(rank = families_order_and_data[i], type = "chronogram_mrca")$tip.label)),
                                error = function(e) paste(print(families_order_and_data[i]))
    )
  }
  names(list_family) <- families_order_and_data
  monotipic_family <- names(unlist(lapply(list_family, function(x) which(length(x) == 1))))
  list_monotipic <- vector(mode = "list", length = length(monotipic_family))
  
  for(i in 1:length(monotipic_family)){
    list_monotipic[[i]]<- tryCatch(fishtree::fishtree_taxonomy(rank = monotipic_family[i])[[1]]$taxonomy[[9]],
                                   error = function(e) paste("not.found", "_", monotipic_family[i], sep = ""))
  } # list of orders which the families belong
  
  orders_to_add <- unique(unlist(list_monotipic[- which(sub("_.*", "", unlist(list_monotipic)) == "not.found")]))
  differences_orders_toadd <- setdiff(rank_order, orders_to_add)
  if(length(differences_orders_toadd) >= 1){
    all_orders_include <- c(differences_orders_toadd, orders_to_add)
  }
  all_orders_include <- unique(c(rank_order, unique(orders_to_add)))
  
  list_order <- vector(mode = "list", length= length(all_orders_include))
  #filtering all species names within orders
  for(i in 1:length(all_orders_include)){
    list_order[[i]]<- tryCatch(paste(print(fishtree::fishtree_phylogeny(rank = all_orders_include[i], type = "chronogram_mrca")$tip.label)),
                               error = function(e) paste(print(all_orders_include[i]))
    )
  }
  names(list_order) <- all_orders_include
  
  # baixando a filogenia com todas as ordens
  
  phylo_order <- filter_rank(order = list_order)  #phylogeny with all order
  phylo_order <- ape::makeNodeLabel(phy = phylo_order) #name nodes for all species
  # phylo_family <- suppressWarnings(filter_rank(order = list_family)) #phylogeny for all family
  
  
  # naming node according to order
  for (i in 1:length(list_order)) {
    # i = 5
    phylo_order<- ape::makeNodeLabel(phylo_order, "u", nodeList = list(Ord_name = list_order[[i]]))
    phylo_order$node.label[which(phylo_order$node.label == "Ord_name")] <- names(list_order)[i]
  }
  
  list_non_monotipic <- list_family[setdiff(names(list_family), monotipic_family)]
  
  # naming node according to families that are not monotipic
  for (i in 1:length(list_non_monotipic)) {
    #i= 225
    phylo_order <- ape::makeNodeLabel(phylo_order, "u", nodeList = list(Fam_name = list_non_monotipic[[i]]))
    phylo_order$node.label[which(phylo_order$node.label == "Fam_name")] <- paste(names(list_non_monotipic)[i])
  }
  
  families_in_tree <- families_order_and_data[which(!is.na(match(families_order_and_data, phylo_order$node.label)) == T)]
  families_monotipic_notfound <- setdiff(monotipic_family, families_in_tree)
  
  for(i in 1:length(families_monotipic_notfound)){
    # i = 73
    
    spp_tmp <-  tryCatch(fishtree::fishtree_taxonomy(rank = families_monotipic_notfound[i])[[1]]$species,
                                                         error = function(e) paste("not.found", "_", families_monotipic_notfound[i], sep = ""))
    spp_tmp <- gsub("\\ ", "_", spp_tmp)
    
    list_family[families_monotipic_notfound[i]] <- list(spp_tmp)
    
  }
  
  phylo_order <- suppressWarnings(filter_rank(order = list_family)) #phylogeny for all family
  phylo_order <- ape::makeNodeLabel(phy = phylo_order)
  
  # repeating the naming procedure to add orders to node labels
  for (i in 1:length(list_order)) {
    # i = 5
    phylo_order<- ape::makeNodeLabel(phylo_order, "u", nodeList = list(Ord_name = list_order[[i]]))
    phylo_order$node.label[which(phylo_order$node.label == "Ord_name")] <- names(list_order)[i]
  }
  
  # ultrametrize the tree
  phylo_order <- phytools::force.ultrametric(phylo_order)
  
  # filtering only for families that was found in the fishtree of life with some available phylogeny
  families_not_found_fishtree <- names(unlist(lapply(
    lapply(list_family, 
           function(x){
             sub("_.*", "", x)
             }
           ), 
    function(y) which(length(y) == 1) & which(y == "not.found")
    )
    )
    )
  
  list_family_tobeaddnames <- list_family[- match(families_not_found_fishtree, names(list_family))]
  
  # repeating the naming procedure to add families to node labels
  for (i in 1:length(list_family_tobeaddnames)) {
    # i = 31
    na_check <- sum(!is.na(match(list_family_tobeaddnames[[i]], phylo_order$tip.label)))
    if(na_check == 1){
      spp_singleton <- unlist(list(list_family_tobeaddnames[[i]][!is.na(match(list_family_tobeaddnames[[i]], phylo_order$tip.label))]))
      spp_singleton_add <- paste(sub("_.*", "", spp_singleton), "_", "singleton", sep = "")
      phylo_order <- phytools::add.species.to.genus(tree = phylo_order, species = spp_singleton_add)
      list_family_tobeaddnames[i] <- list(c(spp_singleton, spp_singleton_add))
    }
    phylo_order <- ape::makeNodeLabel(phylo_order, "u", nodeList = list(Fam_name = list_family_tobeaddnames[[i]]))
    phylo_order$node.label[which(phylo_order$node.label == "Fam_name")] <- paste(names(list_family_tobeaddnames)[i])
    
  }
  
  families_with_no_spp_tree <- names(list_family_tobeaddnames[which(is.na(match(names(list_family_tobeaddnames), phylo_order$node.label)) == TRUE)])
  
  #selecting species that must be added to genus in the tree (sister species)
  spp_data <- 1:length(spp)
  names(spp_data) <- spp
  insert_spp <- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree # all species to be inserted
  phylo_order <- phytools::force.ultrametric(tree = phylo_order)
  
  #species that must be added in the first step of the procedure
  genus_in_tree <- sub("_.*", "", 
                        phylo_order$tip.label)[match(sub("_.*", "", 
                                                         insert_spp), 
                                                     sub("_.*", "", 
                                                         phylo_order$tip.label))][!is.na(sub("_.*", "", 
                                                                                             phylo_order$tip.label)[match(sub("_.*", "", 
                                                                                                                              insert_spp), 
                                                                                                                          sub("_.*", "", 
                                                                                                                              phylo_order$tip.label)
                                                                                             )
                                                                                             ]
                                                         )
                                                         ]
  species_to_genus <- insert_spp[which(is.na(insert_spp[match(sub("_.*", "", insert_spp), genus_in_tree)]) == FALSE)] #genus that must be added
  
  
  #######solving problem 1, first step - add genus
  for(i in 1:length(species_to_genus)){
    #adding species to genus that already exist in the tree
    phylo_order<- phytools::add.species.to.genus(tree = phylo_order, species = species_to_genus[i]) 
  }
  
  ###preparing data to the remaining insertions
  insert_spp2 <- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree #species that must be added after step 1
  
  #second step - adding species that are not genus in the tree, but are in any family
  data_exRound2 <- data[match(insert_spp2, as.character(data$s)),] #data to be submited to round 2 of family search
  rank_family2 <- unique(as.character(data[match(insert_spp2, as.character(data$s)),2]))
  list_spp_step2 <- vector(mode = "list", length= length(rank_family2))
  
  for(i in 1:length(rank_family2)){
    list_spp_step2[[i]]<- tryCatch(paste(ape::extract.clade(phy = phylo_order, node = as.character(rank_family2[i]))$tip.label),
                                   error = function(e) paste("noFamily", as.character(data[which(rank_family2[i] == data$f), 1]), sep= "_"))
  }
  names(list_spp_step2) <- rank_family2
  
  data_exRound3 <- data_exRound2[which(names(which(unlist(lapply(lapply(list_spp_step2, 
                                                                        function(x) which(sub("_.*", "", x) == "noFamily")
  ), 
  function(y) length(y)
  )
  ) > 0)) == data_exRound2$f),] #data to be submited to third round in order search - no species with the same family of these species 
  #in the phylogeny
  
  #species of the same family that species that must be added that are already on the tree
  spp_family <- 1:nrow(data_exRound2)
  names(spp_family) <- data_exRound2$s
  
  spp_with_family <- names(which(unlist(lapply(lapply(list_spp_step2, 
                                                      function(x) which(sub("_.*", "", x) != "noFamily")
  ), 
  function(y) which(length(y) != 0)
  )
  ) > 0))
  
  spp_family_inTree <- suppressWarnings(list_spp_step2[which(names(list_spp_step2) == spp_with_family)])
  
  
  #species to be added in step 2 - species with family representatives
  
  spp_to_add_round2 <- setdiff(data_exRound2$s, data_exRound3$s)
  
  
  ####initializing the insertion of species that present representatives species in family level
  family_name <- names(list_spp_step2)[names(list_spp_step2) == data[which(spp_to_add_round2[1] == data$s), ][, 2]]
  spp_Byfamily_inTree <- as.character(unlist(list_spp_step2[c(which(names(list_spp_step2) == data[which(spp_to_add_round2[1] == data$s), ][, 2]))])) 
  user_option_spp <-  unique(sub("_.*", "", as.character(unlist(spp_Byfamily_inTree))))
  
  if(any(user_option_spp == 1)){
    phylo_order<- phytools::add.species.to.genus(tree = phylo_order, 
                                                 species = paste(sub("_.*", "", as.character(spp_family_inTree)[1])
                                                                 , "toadd", sep= "_"
                                                 )
    )
    position_problem1 <- which(phylo_order$tip.label == paste(sub("_.*", "", as.character(spp_family_inTree)[1]),
                                                              "toadd", sep= "_")
    )
    phylo_test$tip.label[position_problem1]<- spp_to_add_round2 #solving family add when there is only one species of the same family that species to add
  } else{
    while(length(spp_to_add_round2) >= 1){
      family_name <- names(list_spp_step2)[names(list_spp_step2) == data[which(spp_to_add_round2[1] == data$s), ][, 2]]
      spp_Byfamily_inTree <- as.character(unlist(list_spp_step2[c(which(names(list_spp_step2) == data[which(spp_to_add_round2[1] == data$s), ][, 2]))])) 
      user_option_spp <-  unique(sub("_.*", "", as.character(unlist(spp_Byfamily_inTree))))
      local_to_add_spp <- readline(prompt = print_cat(print_cat = user_option_spp, spp = spp_to_add_round2[1], family_name)) #user interactive option to choose species
      
      spp_user_opt <- unlist(strsplit(local_to_add_spp, split = " "))
      
      if(length(spp_user_opt) > 2){
        stop("/n Only two genus can be chosen to insert species")
      }
      if(length(spp_user_opt) < 1){
        stop("/n At least one genus/family can be chosen to insert species")
      }
      
      # if user decide to insert species as politomy between two different genus 
      if(length(spp_user_opt) == 2){
        spp_to_add_tmp <- spp_Byfamily_inTree[match(spp_user_opt, 
                                                    sub("_.*", "", as.character(unlist(spp_Byfamily_inTree)))
        )
        ]
        node_btw_genus <- phytools::fastMRCA(tree = phylo_order, sp1 = spp_to_add_tmp[1], sp2 = spp_to_add_tmp[2])
        phylo_order <- phytools::bind.tip(tree = phylo_order, tip.label = spp_to_add_round2[1], where = node_btw_genus, position = 0)
      }
      
      # if user decided to insert species as politomy of a specific genus or in the family node
      if(length(spp_user_opt) == 1){
        if(length(match(spp_user_opt, c(user_option_spp, family_name))) != 1){
          stop(paste("\n Check the spelling of Genus or Family in:", spp_user_opt, sep = " "))
        }
        if(any(spp_user_opt == family_name)){  # if user decided to insert the species as politomy in family node
          node_family <- which(phylo_order$node.label == family_name)
          phylo_order <- phytools::bind.tip(tree = phylo_order, tip.label = spp_to_add_round2[1], where = node_family, position = 0)
        } else{ 
          genus_nspp <- match(spp_user_opt, names(which(table(sub("_.*", "", as.character(unlist(spp_Byfamily_inTree)))) > 1)))
          if(length(genus_nspp) == 1){
            phylo_order <- ape::makeNodeLabel(phy = phylo_order, method = "u", nodeList = list(MRCA = spp_user_opt))
            position_MRCA <- which(c(phylo_order$tip.label, phylo_order$node.label) == "MRCA")
            size_branch <- phylo_order$edge.length[sapply(position_MRCA, function(x, y) which(y == x), y = phylo_order$edge[, 2])]
            phylo_order <- bind.tip(phylo_order, spp_to_add_round2[1], where = position_MRCA, position = size_branch/2)
          }
          
          if(length(genus_nspp) < 1){ # add as a sister group 
            phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                          species = paste(sub("_.*", "", as.character(spp_user_opt))
                                                                          , "toadd", sep= "_"
                                                          )
            )
            position_problem2 <- which(phylo_order$tip.label == paste(sub("_.*", "", as.character(spp_user_opt)),
                                                                      "toadd", sep= "_")
            )
            phylo_order$tip.label[position_problem2] <- spp_to_add_round2[1]
          }
        }
      }
      
      # running again the check procedure
      
      spp_data <- 1:nrow(data_exRound2)
      names(spp_data) <- data_exRound2$s
      insert_spp <- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree # remaining species to be inserted
      phylo_order<- phytools::force.ultrametric(tree = phylo_order)
      
      # genus checking procedure
      genus_in_tree<- sub("_.*", "", 
                           phylo_order$tip.label)[match(sub("_.*", "", 
                                                            insert_spp), 
                                                        sub("_.*", "", 
                                                            phylo_order$tip.label))][!is.na(sub("_.*", "", 
                                                                                                phylo_order$tip.label)[match(sub("_.*", "", 
                                                                                                                                 insert_spp), 
                                                                                                                             sub("_.*", "", 
                                                                                                                                 phylo_order$tip.label)
                                                                                                )
                                                                                                ]
                                                            )
                                                            ]
      # adding species to genus
      if(length(genus_in_tree) >= 1){
        species_to_genus<- insert_spp[which(is.na(insert_spp[match(sub("_.*", "", insert_spp), genus_in_tree)]) == FALSE)] #genus that must be added
        for(i in 1:length(species_to_genus)){
          #adding species to genus that already exist in the tree
          phylo_order<- phytools::add.species.to.genus(tree = phylo_order, species = species_to_genus[i]) 
        }
        insert_spp<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree # remaining species to be inserted
        phylo_order<- phytools::force.ultrametric(tree = phylo_order)
      }
      
      rank_family2 <- unique(as.character(data[match(insert_spp, as.character(data$s)),2]))
      list_spp_step2 <- vector(mode = "list", length= length(rank_family2))
      for(i in 1:length(rank_family2)){
        list_spp_step2[[i]]<- tryCatch(paste(ape::extract.clade(phy = phylo_order, node = rank_family2[i])$tip.label), 
                                       error = function(e) paste("noFamily", as.character(data[which(rank_family2[i] == data$f), 1]), 
                                                                 sep= "_")
        )
        
      }
      names(list_spp_step2) <- rank_family2
      
      names_family <- unlist(lapply(lapply(list_spp_step2, 
                                           function(x) 
                                             which(sub("_.*", "", x) == "noFamily")
      ), 
      function(y) length(y)
      )
      )
      family_still_add <- names(list_spp_step2)[ifelse(names_family == 1, TRUE, FALSE)] # family that will be added in the next round
      data_exRound3 <- data[match(family_still_add, data$f), ]
      spp_to_add_round2 <- setdiff(insert_spp, data_exRound3$s)
    }
  }
  
  ###### step 3 - add species to orders ######
  has_cichli <- which(data_exRound3$o == "Cichliformes") # There are any Cichliformes in the remaining species
  if(length(has_cichli) >= 1){
    data_exRound3[has_cichli, "o"] <- "Perciformes"
  }
  
  rank_order_Round3 <- rank_order[match(data_exRound3$o, rank_order)]
  families_round3 <- lapply(lapply(rank_order_Round3, function(x){
    fishbase[which(x == fishbase$Order), 10]
  }), function(y) unique(y))
  names(families_round3) <- rank_order_Round3
  
  # initiating insertion in families
  for(i in 1:length(families_round3)){
    i = 1
    user_option_family <-  families_round3[[i]]
    
    local_to_add_spp_family <- readline(prompt = print_cat(print_cat = unlist(user_option_family), 
                                                           spp = data_exRound3$s[i], 
                                                           data_exRound3$o[i])
                                        ) #user interactive option to choose species
    family_user_opt <- unlist(strsplit(local_to_add_spp_family, split = " "))
    
    if(length(family_user_opt) == 1){
      if(family_user_opt == data_exRound3$o[i]){ # insert species in order node
        node_order_pos <- which(phylo_order$node.label == data_exRound3$o[i])
        phytools::bind.tip(tree = phylo_order, 
                           tip.label = data_exRound3$s[i],
                           where = node_order_pos,
                           position = 0)
      } else{
        family_nspp <- length(list_family[[match(family_user_opt, names(list_family))]])
        if(family_nspp > 1){
          position_family <- which(phylo_order$node.label == family_user_opt)
          size_branch_family <- phylo_order$edge.length[sapply(position_family, function(x, y) which(y == x), y = phylo_order$edge[, 2])]
          phylo_order <- phytools::bind.tip(phylo_order, data_exRound3$s[i], where = position_family, position = size_branch_family/2)
        }
        if(family_nspp == 1){
          phylo_order <- phytools::add.species.to.genus(tree = phylo_order, 
                                                        species = paste(sub("_.*", as.character(list_family[[which(names(list_family) == family_user_opt)]][1])),
                                                                        "toadd",
                                                                        sep = "_")
                                                        )
          position_problem_family <- which(phylo_order$tip.label == paste(sub("_.*",
                                                                              "", 
                                                                              as.character(list_family[[which(names(list_family) == family_user_opt)]][1]),
                                                                              "toadd",
                                                                              sep = "_")
                                                                          )
                                           )
          phylo_order$tip.label[position_problem_family] <- data_exRound3$s[i]
        }
      }
    }
  }
  
  #final data proccessing - cutting the phylogenetic tree to obtain only species in data
  data_final<- 1:length(as.character(data$s))
  names(data_final)<- as.character(data$s)
  tree_res<- suppressWarnings(ape::drop.tip(phy = phylo_order, 
                                            tip = treedata_modif(phy = phylo_order, data = data_final)$nc$tree_not_data)
  )
  tree_res #phylogeny with only species on data
}

