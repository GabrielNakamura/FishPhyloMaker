#### internal functions 

# check genus in tree
spp<- as.character(data$s) #especies
spp_data<- 1:length(spp)
names(spp_data)<- spp
insert_spp<- treedata_modif(phy = phylo_order, data = spp_data, warnings = F)$nc$data_not_tree # all species to be inserted
phylo_order<- phytools::force.ultrametric(tree = phylo_order)
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

for(i in 1:length(species_to_genre)){
  #adding species to genus that already exist in the tree
  phylo_order<- phytools::add.species.to.genus(tree = phylo_order, species = species_to_genre[i]) 
}
