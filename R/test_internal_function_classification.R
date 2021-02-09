level <- "s"
data_to_insert <- insert_spp2
data <- taxon_data
data <- edit(taxon_data)

function(data, data_to_insert, level= "s"){
  
  l1 <- colnames(data)[match(level, colnames(data))]
  pmatch(level, colnames(data))
  data_exRound2 <- data[match(data_to_insert, as.character(data[, l1])),] #data to be submited to round 2 of family search
  rank_family2 <- as.character(data[match(data_to_insert, as.character(data[, l1])),2])
  list_spp_step2<- vector(mode = "list", length= length(rank_family2))
  
  for(i in 1:length(rank_family2)){
    list_spp_step2[[i]]<- tryCatch(paste(ape::extract.clade(phy = phylo_order, node = as.character(rank_family2[i]))$tip.label),
                                   error = function(e) paste("noFamily", as.character(data[which(rank_family2[i] == data[, pmatch(level, colnames(data)) + 1]), 1]), sep= "_"))
  }
  
  
  data_exRound3<- data_exRound2[which(unlist(lapply(lapply(list_spp_step2, 
                                                           function(x) which(sub("_.*", "", x) == "noFamily")
  ), 
  function(y) length(y)
  )
  ) > 0),] #data to be submited to third round in order search - no species with the same family of these species 
  #in the phylogeny
  
  data_exRoundFamily <- data[unique(unlist(lapply(as.character(data_exRound2[, pmatch(level, colnames(data)) + 1]), 
                                                 function(x) which(x == as.character(data[, pmatch(level, colnames(data)) + 1])
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
  user_option_spp<-  unique(sub("_.*", "", as.character(spp_family_inTree))) # genus in tree
  
  ####initializing the insertion of species that present representatives species in family level
  if(user_option_spp == 1){
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
      
      # if user decide to insert species as politomy between two different genus 
      if(length(local_to_add_spp) == 2){
        local_to_add_spp_btw <- unlist(strsplit(local_to_add_spp, split = " "))
        spp_to_add_tmp <- spp_family_inTree[match(local_to_add_spp_btw, sub("_.*", "", as.character(spp_family_inTree)))]
        node_btw_genus <- phytools::fastMRCA(tree = phylo_order, sp1 = spp_to_add_tmp[1], sp2 = spp_to_add_tmp[2])
        phylo_order <- phytools::bind.tip(tree = phylo_order, tip.label = spp_to_add_round2[1], where = node_btw_genus, position = 0)
      }
      
      # if user decided to insert species as politomy of a specific genus
      if(length(local_to_add_spp) == 1){
        phylo_order<- phytools::add.species.to.genus(tree = phylo_order, 
                                                     species = paste(sub("_.*", "", as.character(local_to_add_spp))
                                                                     , "toadd", sep= "_"
                                                     )
        )
        position_problem2<- which(phylo_order$tip.label == paste(sub("_.*", "", as.character(local_to_add_spp)),
                                                                 "toadd", sep= "_")
        )
        phylo_order$tip.label[position_problem2]<- spp_to_add_round2[1] 
        
      }
      
      # if user decided to insert the species as politomy in family node
      if(local_to_add_spp == "family_node"){
        node_family <- which(phylo_order$node.label ==  data[match(spp_to_add_round2[1], data$s), "f"])
        phylo_order <- phytools::bind.tip(tree = phylo_order, tip.label = spp_to_add_round2[1], where = node_family, position = 0)
      }
      
      #running again the check procedure
      rank_family2 <- as.character(data[match(data_to_insert, as.character(data$s)),2])
      data_to_insert<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree #species that must be added after step 1
      data_exRound2<- data[match(data_to_insert, as.character(data$s)),]
      list_spp_step2<- vector(mode = "list", length= length(unique(rank_family2)))
      for(i in 1:length(unique(rank_family2))){
        #i= 1
        list_spp_step2[[i]]<- tryCatch(paste(ape::extract.clade(phy = phylo_order, node = as.character(unique(rank_family2)[i]))$tip.label), 
                                       error = function(e) paste("noFamily", as.character(data[which(unique(rank_family2)[i] == data$f), 1]), 
                                                                 sep= "_")
        )
        
      }
      
      data_exRound3<- data_exRound2[1:unlist(lapply(lapply(list_spp_step2, 
                                                           function(x) 
                                                             which(sub("_.*", "", x) == "noFamily")
      ), 
      function(y) length(y)
      )
      ),
      ] #data to be submited to third round in order search - no species with the same family of these species 
      #in the phylogeny
      spp_to_add_round2<- setdiff(data_exRound2$s, data_exRound3$s)
    }
  }
}
