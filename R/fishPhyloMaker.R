#' Making fish phylogeny for a local pool of species 
#'
#' @param data A character vector containing the names of species or a dataframe generated from \code{tab_function}
#'
#' @return A newick object describing the phylogenetic relationships among the species provided in \code{data}
#' 
#' @export
#'
#' @examples
#' 
phyloMatch<- function(data){
  
  #organizing taxonomic levels
  rank_order <- as.character(unique(data$o)) #ordens
  rank_family <- as.character(unique(data$f)) #familias
  spp <- as.character(data$s) #especies
  list_ordem <- vector(mode = "list", length= length(rank_order))
  list_family <- vector(mode = "list", length= length(rank_family))
  
  #filtering all species names within orders presented in the original phylogeny
  for(i in 1:length(rank_order)){
    list_ordem[[i]] <- tryCatch(paste(print(fishtree::fishtree_phylogeny(rank = rank_order[i], type = "chronogram_mrca")$tip.label)),
                               error = function(e) paste(print(rank_order[i]))
    )
  }
  
  #filtering for family
  for(i in 1:length(rank_family)){
    list_family[[i]] <- tryCatch(paste(print(fishtree::fishtree_phylogeny(rank = rank_family[i], type = "chronogram_mrca")$tip.label)),
                                error = function(e) paste(print(rank_family[i]))
    )
  }
  
  ##downloading phylogeny from all orders in data
  filter_rank <- function(ordem){
    if(length(which(sub("_.*", "", unlist(ordem)) == "not.found")) >= 1){
      
      phy_ord<- fishtree::fishtree_phylogeny(unlist(ordem)[-which(sub("_.*", "", unlist(ordem)) == "not.found")])
      
    } else{
      phy_ord<- fishtree::fishtree_phylogeny(unlist(ordem))
    }
    phy_ord
  }
  
  phylo_order <- filter_rank(ordem = list_ordem)  #phylogeny with all order
  phylo_order <- ape::makeNodeLabel(phy = phylo_order) #name nodes for all species
  phylo_family <- suppressWarnings(filter_rank(ordem = list_family)) #phylogeny for all family
  
  
  #naming node according to order
  for( i in 1:length(list_ordem)){
    temp <- list_ordem[[i]]
    phylo_temp <- ape::drop.tip(phy = phylo_order,  setdiff(phylo_order$tip.label, temp))
    node_ordem <- phylo_temp$node.label[1]
    phylo_order$node.label[which(phylo_order$node.label == node_ordem)] <- paste(rank_order[i])
  }
  
  #naming node family in phylo order
  for( i in 1:length(list_family)){
    temp <- list_family[[i]]
    phylo_temp <- suppressWarnings(ape::drop.tip(phy = phylo_order,  setdiff(phylo_order$tip.label, temp)))
    node_ordem <- phylo_temp$node.label[1]
    phylo_order$node.label[which(phylo_order$node.label == node_ordem)] <- paste(rank_family[i])
  }
  
  #selecting species that must be added to genus in the tree (sister species)
  spp_data<- 1:length(spp)
  names(spp_data)<- spp
  insert_spp<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree #all species to be inserted
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
  species_to_genre<- insert_spp[which(is.na(insert_spp[match(sub("_.*", "", insert_spp), genres_in_tree)]) == FALSE)]
  #genre that must be added
  
  
  #######solving problem 1, first step - add genus
  for(i in 1:length(species_to_genre)){
    #adding species to genus that already exist in the tree
    phylo_order<- phytools::add.species.to.genus(tree = phylo_order, species = species_to_genre[i]) 
  }
  
  ###preparing data to the remaining insertions
  insert_spp2<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree #species that must be added after step 1
  
  if(length(insert_spp2) == 0){
    data_final<- 1:length(as.character(data$s))
    names(data_final)<- as.character(data$s)
    tree_res<- suppressWarnings(ape::drop.tip(phy = phylo_order, 
                                              tip = treedata_modif(phy = phylo_order, data = data_final)$nc$tree_not_data)
    )
    return(tree_res)
  } else{
    #second step - adding species that are not genus in the tree, but are in any family 
    data_exRound2<- data[match(insert_spp2, as.character(data$s)),] #data to be submited to round 2 of family search
    rank_family2<- as.character(data[match(insert_spp2, as.character(data$s)),2])
    list_spp_step2<- vector(mode = "list", length= length(rank_family2))
    
    for(i in 1:length(rank_family2)){
      #i= 1
      list_spp_step2[[i]]<- tryCatch(paste(ape::extract.clade(phy = phylo_order, node = as.character(rank_family2[i]))$tip.label),
                                     error = function(e) paste("noFamily", as.character(data[which(rank_family2[i] == data$f), 1]), sep= "_"))
      
    }
    
    
    data_exRound3<- data_exRound2[which(unlist(lapply(lapply(list_spp_step2, 
                                                             function(x) which(sub("_.*", "", x) == "noFamily")
    ), 
    function(y) length(y)
    )
    ) > 0),] #data to be submited to third round in order search - no species with the same family of these species 
    #in the phylogeny
    
    data_exRoundFamily<- data[unique(unlist(lapply(as.character(data_exRound2$f), 
                                                   function(x) which(x == as.character(data$f)
                                                   )
    )
    )
    )
    , ]
    
    
    #species of the same family that species that must be added that are already on the tree
    spp_family<- 1:nrow(data_exRoundFamily)
    names(spp_family)<- data_exRoundFamily$s
    spp_family_inTree<- as.character(data_exRoundFamily$s[-match(as.character(suppressWarnings(treedata_modif(phy = phylo_order, spp_family)$nc$data_not_tree)),
                                                                 as.character(data_exRoundFamily$s)
    )
    ]
    ) 
    #species to be added in step 2 - species with family representatives
    spp_to_add_round2<- setdiff(data_exRound2$s, data_exRound3$s)
    user_option_spp<- spp_family_inTree
    
    #### spp family in supertree but not on data
    family_inSuperTree <- setdiff(unique(data[match(spp_to_add_round2, data$s), "f"]), 
                                  unique(data[match(spp_family_inTree, data$s), "f"]))
    
    spp_family_inSuperTree <- data[match(family_inSuperTree, data$f), "s"]
    
    
    ### updating ssp to add in family round removing spp with family in supertree
    spp_to_add_round2 <- spp_to_add_round2[-match(spp_family_inSuperTree, spp_to_add_round2)]
    
    ####initializing the insertion of species that present representatives species in family level
    if(length(spp_family_inTree) == 1){
      phylo_order<- phytools::add.species.to.genus(tree = phylo_order, 
                                                   species = paste(sub("_.*", "", as.character(spp_family_inTree)[1])
                                                                   , "toadd", sep= "_"
                                                   )
      )
      position_problem1<- which(phylo_order$tip.label == paste(sub("_.*", "", as.character(spp_family_inTree)[1]),
                                                               "toadd", sep= "_")
      )
      phylo_test$tip.label[position_problem1]<- spp_to_add_round2 #solving family add when there is only one species of the same family that species to add
    } else{
      while(length(spp_to_add_round2) >= 1){
        local_to_add_spp<- readline(prompt = print_cat(print_cat = user_option_spp, spp = spp_to_add_round2[1])) #user interactive option to choose species
        phylo_order<- phytools::add.species.to.genus(tree = phylo_order, 
                                                     species = paste(sub("_.*", "", as.character(local_to_add_spp))
                                                                     , "toadd", sep= "_"
                                                     )
        )
        position_problem2<- which(phylo_order$tip.label == paste(sub("_.*", "", as.character(local_to_add_spp)),
                                                                 "toadd", sep= "_")
        )
        phylo_order$tip.label[position_problem2]<- spp_to_add_round2[1] #solving problem 2 and 3 to insert species in family and/or genre
        
        #running again the check procedure
        insert_spp2<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree #species that must be added after step 1
        rank_family2<- as.character(data[match(insert_spp2, as.character(data$s)),2])
        data_exRound2<- data[match(insert_spp2, as.character(data$s)),]
        list_spp_step2<- vector(mode = "list", length= length(rank_family2))
        for(i in 1:length(rank_family2)){
          #i= 1
          list_spp_step2[[i]]<- tryCatch(paste(ape::extract.clade(phy = phylo_order, node = as.character(rank_family2)[i])$tip.label), 
                                         error = function(e) paste("noFamily", as.character(data[which(rank_family2[i] == data$f), 1]), 
                                                                   sep= "_")
          )
          
        }
        
        data_exRound3<- data_exRound2[which(unlist(lapply(lapply(list_spp_step2, 
                                                                 function(x) which(sub("_.*", "", x) == "noFamily")
        ), 
        function(y) length(y)
        )
        ) > 0),] #data to be submited to third round in order search - no species with the same family of these species 
        #in the phylogeny
        
        spp_to_add_round2<- setdiff(data_exRound2$s, data_exRound3$s)
        family_inSuperTree <- setdiff(unique(data[match(spp_to_add_round2, data$s), "f"]), 
                                      unique(data[match(spp_family_inTree, data$s), "f"]))
        spp_family_inSuperTree <- data[match(family_inSuperTree, data$f), "s"]
        spp_to_add_round2 <- spp_to_add_round2[-match(spp_family_inSuperTree, spp_to_add_round2)]
        
        data_exRoundFamily<- data[unique(unlist(lapply(as.character(data_exRound2$f), 
                                                       function(x) which(x == as.character(data$f)
                                                       )
        )
        )
        )
        , ]
        
        
        spp_family<- 1:nrow(data_exRoundFamily)
        names(spp_family)<- data_exRoundFamily$s
        spp_family_inTree<- as.character(data_exRoundFamily$s[-match(as.character(suppressWarnings(treedata_modif(phy = phylo_order, spp_family)$nc$data_not_tree)),
                                                                     as.character(data_exRoundFamily$s)
        )
        ]
        ) 
        
        #species to be added in step 2 - species with family representatives
        user_option_spp<- spp_family_inTree
        spp_to_add_round2<- setdiff(data_exRound2$s, data_exRound3$s)
        spp_to_add_round2 <- spp_to_add_round2[-match(spp_family_inSuperTree, spp_to_add_round2)]
      }
    }
    
    ### adding species to family level but are not on spp pool
    if(length(family_inSuperTree) == 1){
      sample_family_inSuperTree <- sample(tryCatch(paste(ape::extract.clade(phy = phylo_order, node = family_inSuperTree)$tip.label), 
                                                   error = function(e) paste("noFamily", as.character(data[which(family_inSuperTree == data$f), 1]), 
                                                                             sep= "_")
      ), 1)
      spp_to_add_Family <- paste(sub("_.*", "", as.character(sample_family_inSuperTree)[1]), "toadd", sep= "_")
      phylo_order <- phytools::add.species.to.genus(tree = phylo_order, species = spp_to_add_Family)
      
      position_problem_family <- which(phylo_order$tip.label == spp_to_add_Family)
      phylo_order$tip.label[position_problem_family]<- spp_family_inSuperTree #solving family add when there is only one species of the same family that species to add
    } 
    
    # for two species with spp representants in superTree at family level 
    if(length(family_inSuperTree) == 2){
      for(i in 1:length(spp_family_inSuperTree)){
        sample_family_inSuperTree <- sample(tryCatch(paste(ape::extract.clade(phy = phylo_order, node = family_inSuperTree)$tip.label), 
                                                     error = function(e) paste("noFamily", as.character(data[which(family_inSuperTree == data$f), 1]), 
                                                                               sep= "_")
        ), 1)
        spp_to_add_Family <- paste(sub("_.*", "", as.character(sample_family_inSuperTree)[1]), "toadd", sep= "_")
        phylo_order <- phytools::add.species.to.genus(tree = phylo_order, species = spp_to_add_Family)
        position_problem_family <- which(phylo_order$tip.label == spp_to_add_Family)
        phylo_order$tip.label[position_problem_family]<- spp_family_inSuperTree[i] #solving family add when there is only one species of the same family that species to add
      }
    }
    
    # for more than two species with representatives in the 
    if(length(family_inSuperTree) > 3){
      sample_family_inSuperTree <- numeric(length = 2)
      for(i in 1:2){
        sample_family_inSuperTree[i] <- sample(tryCatch(paste(ape::extract.clade(phy = phylo_order, node = family_inSuperTree)$tip.label), 
                                                        error = function(e) paste("noFamily", as.character(data[which(family_inSuperTree == data$f), 1]), 
                                                                                  sep= "_")
        ), 1)
        spp_to_add_Family <- paste(sub("_.*", "", as.character(sample_family_inSuperTree)[1]), "toadd", sep= "_")
        phylo_order <- phytools::add.species.to.genus(tree = phylo_order, species = spp_to_add_Family)
        position_problem_family <- which(phylo_order$tip.label == spp_to_add_Family)
        phylo_order$tip.label[position_problem_family]<- spp_family_inSuperTree[i] #solving family add when there is only one species of the same family that species to add  
        
      }
      user_option_spp_familySupTree <- phylo_order$tip.label[which(phylo_order$tip.label == spp_family_inSuperTree)]
      while(spp_family_inSuperTree >= 1){
        spp_to_add_round2_family<- setdiff(data_exRound2$s, data_exRound3$s)
        local_to_add_spp_family <- readline(prompt = print_cat(print_cat = user_option_spp_familySupTree, spp = spp_to_add_round2_family[1])) #user interactive option to choose species
        spp_to_add_Family <- paste(sub("_.*", "", as.character(spp_to_add_round2_family)[1]), "toadd", sep= "_")
        phylo_order <- phytools::add.species.to.genus(tree = phylo_order, species = spp_to_add_Family)
        position_problem_family <- which(phylo_order$tip.label == spp_to_add_Family)
        phylo_order$tip.label[position_problem_family]<- spp_family_inSuperTree[1] #solving family add when there is only one species of the same family that species to add
        spp_family_inSuperTree[-1]
      }
    }
    ######step 3 - add species to orders######
    species_order_inTree<- match(data$s, 
                                 ape::extract.clade(phy = phylo_order, node = unique(as.character(data_exRound3$o)))$tip.label)[which(!is.na(match(data$s, 
                                                                                                                                                   ape::extract.clade(phy = phylo_order, node = unique(as.character(data_exRound3$o)))$tip.label)) == T)] #species that are already on the tree
    spp_orderTree<- phylo_order$tip.label[species_order_inTree] #species names from orders of species that must be added
    spp_to_add_round3<- as.character(data_exRound3$s) #species that must be added
    if(dim(data_exRound3)[1] >= 1){
      if(sum(species_order_inTree) <= 1){ #is ther any species of this order already inserted in the phylogenetic tree?
        if(dim(data_exRound3)[1] <= 2){
          for(i in 1:dim(data_exRound3)[1]){
            #i= 1
            phylo_order<- phytools::bind.tip(tree = phylo_order, tip.label = as.character(data_exRound3$s)[i], 
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
              phylo_order<- phytools::bind.tip(tree = phylo_order, tip.label = as.character(data_exRound3$s)[k], 
                                               where = which(phylo_order$node.label == as.character(data_exRound3$o[k])) + 1)
            }
          }
        }
      } else{
        for(l in 1:length(spp_to_add_round3)){
          #l= 1
          local_to_add_spp<- readline(prompt = print_cat(print_cat = species_order_inTree, spp = spp_to_add_round3[l]), family = data_exRound3$f[l]) #user interactive option to choose species
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
    return(tree_res) #phylogeny with only species on data
  }
}
