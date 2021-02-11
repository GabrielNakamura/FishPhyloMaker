data_all <- load(here::here("data", "neotropical_comm.rda"))
data_comm <- neotropical_comm[, -c(1, 2)]
source(here::here("R", "tab_function.R"))
source(here::here("R", "internal_user_opt_printCat.R"))
source(here::here("R", "internal_user_opt_printCat2.R"))
source(here::here("R", "internal_treedata_modif.R"))
taxon_data <- tab_function(data_comm)
Loricariidae
Siluriformes
data_process <- taxon_data[c(1, 2, 6, 29, 58), ]
data_process <- rbind(data_process, c("Curculionichthys_inesperado", "Loricariidae", "Siluriformes") )
data_process <- rbind(data_process, c("Dinkiwinki_dipsii", "Teletubiidae", "Cichliformes"))
data_process <- rbind(data_process, c("Peixe_loricariaentregeneros", "Loricariidae", "Siluriformes"))
data_process <- rbind(data_process, c("Peixo_basefamilia", "Loricariidae", "Siluriformes"))
data <- data_process
phyloMatch<- function(data){
  
  #organizing taxonomic levels
  rank_order <- as.character(unique(data$o)) #ordens
  rank_family <- as.character(unique(data$f)) #familias
  spp <- as.character(data$s) #especies
  list_ordem <- vector(mode = "list", length= length(rank_order))
  list_family <- vector(mode = "list", length= length(rank_family))
  
  # list of families within genus presented in data
  cichliformes_ord <- which(rank_order == "Cichliformes")
  if(length(cichliformes_ord) == 1){
    rank_order[cichliformes_ord] <- "Perciformes"
  }
  all_families <- unique(unlist(lapply(rank_order, function(x){
    fishbase[which(x == fishbase$Order), 10]
  })))
  families_in_orders <- suppressWarnings(all_families[which(unique(data$f) != all_families)])
  families_order_and_data <- c(rank_family, families_in_orders)
  
  #filtering all species names within orders
  for(i in 1:length(rank_order)){
    list_ordem[[i]]<- tryCatch(paste(print(fishtree::fishtree_phylogeny(rank = rank_order[i], type = "chronogram_mrca")$tip.label)),
                               error = function(e) paste(print(rank_order[i]))
    )
  }
  
  #filtering for family
  for(i in 1:length(families_order_and_data)){
    list_family[[i]]<- tryCatch(paste(print(fishtree::fishtree_phylogeny(rank = families_order_and_data[i], type = "chronogram_mrca")$tip.label)),
                                error = function(e) paste(print(families_order_and_data[i]))
    )
  }
  
  ##downloading phylogeny from all orders in data
  filter_rank<- function(ordem){
    if(length(which(sub("_.*", "", unlist(ordem)) == "not.found")) >= 1){
      
      phy_ord<- fishtree::fishtree_phylogeny(unlist(ordem)[-which(sub("_.*", "", unlist(ordem)) == "not.found")])
      
    } else{
      phy_ord<- fishtree::fishtree_phylogeny(unlist(ordem))
    }
    phy_ord
  }
  
  phylo_order<- filter_rank(ordem = list_ordem)  #phylogeny with all order
  phylo_order<- ape::makeNodeLabel(phy = phylo_order) #name nodes for all species
  phylo_family<- suppressWarnings(filter_rank(ordem = list_family)) #phylogeny for all family
  
  
  #naming node according to order
  for( i in 1:length(list_ordem)){
    temp<- list_ordem[[i]]
    phylo_temp<- ape::drop.tip(phy = phylo_order,  setdiff(phylo_order$tip.label, temp))
    node_ordem<- phylo_temp$node.label[1]
    phylo_order$node.label[which(phylo_order$node.label == node_ordem)]<- paste(rank_order[i])
  }
  
  #naming node family in phylo order
  for( i in 1:length(list_family)){
    temp<- list_family[[i]]
    phylo_temp<- suppressWarnings(ape::drop.tip(phy = phylo_order,  setdiff(phylo_order$tip.label, temp)))
    node_ordem<- phylo_temp$node.label[1]
    phylo_order$node.label[which(phylo_order$node.label == node_ordem)]<- paste(families_order_and_data[i])
  }
  
  #selecting species that must be added to genus in the tree (sister species)
  spp_data<- 1:length(spp)
  names(spp_data)<- spp
  insert_spp<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree # all species to be inserted
  phylo_order<- phytools::force.ultrametric(tree = phylo_order)
  
  #species that must be added in the first step of the procedure
  genres_in_tree<- sub("_.*", "", 
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
  species_to_genre<- insert_spp[which(is.na(insert_spp[match(sub("_.*", "", insert_spp), genres_in_tree)]) == FALSE)] #genus that must be added
  
  
  #######solving problem 1, first step - add genus
  for(i in 1:length(species_to_genre)){
    #adding species to genus that already exist in the tree
    phylo_order<- phytools::add.species.to.genus(tree = phylo_order, species = species_to_genre[i]) 
  }
  
  ###preparing data to the remaining insertions
  insert_spp2<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree #species that must be added after step 1
  
  #second step - adding species that are not genus in the tree, but are in any family
  data_exRound2<- data[match(insert_spp2, as.character(data$s)),] #data to be submited to round 2 of family search
  rank_family2 <- unique(as.character(data[match(insert_spp2, as.character(data$s)),2]))
  list_spp_step2<- vector(mode = "list", length= length(rank_family2))
  
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
  spp_family<- 1:nrow(data_exRound2)
  names(spp_family)<- data_exRound2$s
  
  spp_with_family <- names(which(unlist(lapply(lapply(list_spp_step2, 
                                   function(x) which(sub("_.*", "", x) != "noFamily")
  ), 
  function(y) which(length(y) != 0)
  )
  ) > 0))
  
  spp_family_inTree <- list_spp_step2[which(names(list_spp_step2) == spp_with_family)]
  
   
  #species to be added in step 2 - species with family representatives
  
  spp_to_add_round2 <- setdiff(data_exRound2$s, data_exRound3$s)
  
  
  ####initializing the insertion of species that present representatives species in family level
  family_name <- names(list_spp_step2)[names(list_spp_step2) == data[which(spp_to_add_round2[1] == data$s), ][, 2]]
  spp_Byfamily_inTree <- as.character(unlist(list_spp_step2[c(which(names(list_spp_step2) == data[which(spp_to_add_round2[1] == data$s), ][, 2]))])) 
  user_option_spp <-  unique(sub("_.*", "", as.character(unlist(spp_Byfamily_inTree))))
  
  if(user_option_spp == 1){
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
        if(any(spp_user_opt == family_name)){  # if user decided to insert the species as politomy in family node
          node_family <- which(phylo_order$node.label == family_name)
          phylo_order <- phytools::bind.tip(tree = phylo_order, tip.label = spp_to_add_round2[1], where = node_family, position = 0)
        } else{ # if user decided to insert species as politomy of a specific genus
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
      
      # running again the check procedure
      
      
      spp_data<- 1:nrow(data_exRound2)
      names(spp_data)<- data_exRound2$s
      insert_spp<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree # remaining species to be inserted
      phylo_order<- phytools::force.ultrametric(tree = phylo_order)
      
      # genus checking procedure
      genres_in_tree<- sub("_.*", "", 
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
      if(length(genres_in_tree) >= 1){
        species_to_genre<- insert_spp[which(is.na(insert_spp[match(sub("_.*", "", insert_spp), genres_in_tree)]) == FALSE)] #genus that must be added
        for(i in 1:length(species_to_genre)){
          #adding species to genus that already exist in the tree
          phylo_order<- phytools::add.species.to.genus(tree = phylo_order, species = species_to_genre[i]) 
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
      spp_to_add_round2<- setdiff(insert_spp, data_exRound3$s)
    }
  }
  
  ######step 3 - add species to orders######
  #naming node according to order
  for( i in 1:length(list_ordem)){
    temp<- list_ordem[[i]]
    phylo_temp<- ape::drop.tip(phy = phylo_order,  setdiff(phylo_order$tip.label, temp))
    node_ordem<- phylo_temp$node.label[1]
    phylo_order$node.label[which(phylo_order$node.label == node_ordem)]<- paste(rank_order[i])
  }
  
  for(i in 1:length(rank_order)){
    i = 4
    list_ordem[[i]]<- tryCatch(paste(print(fishtree::fishtree_phylogeny(rank = rank_order[i], type = "chronogram_mrca")$tip.label)),
                               error = function(e) paste(print(rank_order[i]))
    )
  }
  test_perciformes <- fishtree::fishtree_phylogeny(rank = rank_order[i], type = "chronogram_mrca")
  test_perciformes
  test_perciformes <- ape::makeNodeLabel(phy = test_perciformes)
  names_perciformes <- test_perciformes$tip.label
  
  
  node_ordem<- phylo_temp$node.label[1]
  for( i in 1:length(list_family)){
    temp<- list_family[[i]]
    phylo_temp<- suppressWarnings(ape::drop.tip(phy = phylo_order,  setdiff(phylo_order$tip.label, temp)))
    node_ordem<- phylo_temp$node.label[1]
    phylo_order$node.label[which(phylo_order$node.label == node_ordem)]<- paste(families_order_and_data[i])
  }
  
  #naming node family in phylo order
  
  check_cichli <- which(data_exRound3$o == "Cichliformes")
  
  if(length(check_cichli) >= 1){
    data_exRound3[which(data_exRound3$o == "Cichliformes"), "o"] <- "Perciformes"
  }
  
  for(i in 1:length(data_exRound3$o)){
    i = 1
    list_families_in_order <- ape::extract.clade(phy = phylo_order, 
                       node = unique(as.character(data_exRound3$o[i])))$node.label
    list_families_in_order[which(!is.na(match(families_order_and_data, list_families_in_order)) == TRUE)]
    length(list_families_in_order)
  }
  length(list_families_in_order)
  length(families_in_orders)
  match(families_order_and_data, list_families_in_order)
  
  
  species_order_inTree <- match(data$s, 
                              ape::extract.clade(phy = phylo_order, 
                                                 node = unique(as.character(data_exRound3$o)))$tip.label
                              )[which(!is.na(match(data$s, 
                                                   ape::extract.clade(phy = phylo_order, 
                                                                      node = unique(as.character(data_exRound3$o)))$tip.label)
                                             ) == T)] #species that are already on the tree
  spp_orderTree<- phylo_order$tip.label[species_order_inTree] #species names from orders of species that must be added
  spp_to_add_round3<- as.character(data_exRound3$s) #species that must be added
  if(dim(data_exRound3)[1] >= 1){
    if(sum(species_order_inTree) <= 1){ # is there any species of this order already inserted in the phylogenetic tree?
      if(dim(data_exRound3)[1] <= 2){
        for(i in 1:dim(data_exRound3)[1]){
          #i= 1
          phylo_order<- bind.tip(tree = phylo_order, tip.label = as.character(data_exRound3$s)[i], 
                                 where = which(phylo_order$node.label == as.character(data_exRound3$o[i])) + 1)
        } 
      } else{
        #if there is only one order that species must be added
        user_option_spp2<- print_cat2(spp_to_add_round3)
        if(ape::is.phylo(user_option)){
          phylo_order<- ape::bind.tree(x = phylo_order, y = user_option_spp2, 
                                       where = which(phylo_order$node.label == as.character(data_exRound3$o[1])) + 1) #adding the multiple species according to newick file
        }
        if(user_option_spp2 == "politomy"){
          for(k in 1:length(spp_to_add_round3)){
            phylo_order<- bind.tip(tree = phylo_order, tip.label = as.character(data_exRound3$s)[k], 
                                   where = which(phylo_order$node.label == as.character(data_exRound3$o[k])) + 1)
          }
        }
      }
    } else{
      for(l in 1:length(spp_to_add_round3)){
        #l= 1
        local_to_add_spp<- readline(prompt = print_cat(print_cat = species_order_inTree, spp = spp_to_add_round3[l])) #user interactive option to choose species
        phylo_order<- phytools::add.species.to.genus(tree = phylo_order, 
                                                     species = paste(sub("_.*", "", as.character(local_to_add_spp))
                                                                     , "toadd", sep= "_"
                                                     )
        )
        position_problem3<- which(phylo_order$tip.label == paste(sub("_.*", "", as.character(local_to_add_spp)),
                                                                 "toadd", sep= "_")
        )
        phylo_order$tip.label[position_problem3]<- spp_to_add_round3[l] #solving problem 2 and 3 to insert species in family and/or genre
      }
    }
  } else{
    tree_res<- ape::drop.tip(phy = phylo_order, tip = treedata_modif(phy = phylo_order, data = data$s)$nc$data_not_tree)
  }
  
  #final data proccessing - cutting the phylogenetic tree to obtain only species in data
  data_final<- 1:length(as.character(data$s))
  names(data_final)<- as.character(data$s)
  tree_res<- suppressWarnings(ape::drop.tip(phy = phylo_order, 
                                            tip = treedata_modif(phy = phylo_order, data = data_final)$nc$tree_not_data)
  )
  tree_res #phylogeny with only species on data
}
