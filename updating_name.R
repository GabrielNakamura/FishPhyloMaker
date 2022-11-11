# Data and packages

library(ape)
library(devtools)
library(phytools)
load_all()
data <- spp_df <- read.table('/Users/gabriel.nakamuradesouza/Library/CloudStorage/OneDrive-Personal/Manuscritos/Darwinian_shortfalls/data/taxa_table.txt', header = TRUE)
data_insertion <- data
data_insertion$insertion <- NA


# download whole phylogeny 

tree <- fishtree::fishtree_phylogeny()
tree <- ape::makeNodeLabel(phy = tree)
tax <- fishtree::fishtree_taxonomy()
families <- tax[which(tax$rank == "family"), ]
orders <-  tax[which(tax$rank == "order"), ]
spp <- lapply(families$name, function(x) fishtree::fishtree_taxonomy(ranks = x)[[1]]$sampled_species)
names(spp) <- families$name
spp <- lapply(spp, function(x) gsub(" ", "_", x))


# naming families with more than one species phylogeny
not_samp_families <- names(unlist(lapply(spp, function(x) which(length(x) == 0)))) # fdmilies not sampled in phylogeny
monotip_families <-  names(unlist(lapply(spp, function(x) which(length(x) == 1)))) # monotipic families in the tree
one_spp <- c(not_samp_families, monotip_families)
family_names_tree <- spp[-match(one_spp, names(spp))]
position_family <- 
  lapply(family_names_tree, function(x){
    node_anc <- ape::getMRCA(phy = tree, tip = x)
  })
for(i in 1:length(position_family)){ # naming tree nodes
  tree$node.label[position_family[[i]] - ape::Ntip(tree)] <- names(position_family)[i]
}

# adding fake species to name monotipic families in backbone tree
tree <- phytools::force.ultrametric(tree)
tree_update <- tree
monotipic_list <- spp[match(monotip_families, names(spp))]

for(i in 1:length(monotip_families)){ 
  # i = 1
  spp_monotipic <- spp[match(monotip_families[i], names(spp))][[1]]
  
  spp_singleton <- tree$tip.label[match(spp_monotipic, 
                                        tree$tip.label)]
  spp_singleton_add <- paste(sub("_.*", "", spp_singleton), 
                             "_", "singleton", sep = "")
  tree_update <- phytools::add.species.to.genus(tree = tree_update, 
                                            species = spp_singleton_add)
  monotipic_list[[i]] <- c(spp_singleton, spp_singleton_add)
}

# naming monotipic families

position_monotipic <- 
  lapply(monotipic_list, function(x){
    node_anc <- ape::getMRCA(phy = tree_update, tip = x)
  })
for(i in 1:length(position_monotipic)){ # naming tree nodes
  tree_update$node.label[position_monotipic[[i]] - ape::Ntip(tree_update)] <- names(position_monotipic)[i]
}

# species already in tree
spp_in_tree <- data[data$s %in% tree_update$tip.label, "s"]

data_insertion[match(spp_in_tree, data_insertion$s), "insertion"] <- "present in tree"


# congeneric grafting ----------------------------------------------------

tree_update_genus <- tree_update
pb_congeneric <- progress::progress_bar$new(format = "Adding congeneric species [:bar] :percent", 
                                            total = length(insert_genus), clear = FALSE, 
                                            width = 60, current = "<", incomplete = ">", 
                                            complete = ">")

data_genus <- data[-which(data$s %in% spp_in_tree == TRUE), ]
genus_data <- sub("_.*", "", data_genus$s)
genus_tree <- sub("_.*", "", tree_update_genus$tip.label)
genus_in_tree <- genus_data[genus_data %in% genus_tree]
insert_genus <- data_genus[genus_data %in% genus_tree, "s"]
pb_congeneric <- progress::progress_bar$new(format = "Adding congeneric species [:bar] :percent", 
                                            total = length(insert_genus), clear = FALSE, 
                                            width = 60, current = "<", incomplete = ">", 
                                            complete = ">")
for(i in 1:length(insert_genus)){
  tree_update_genus <- phytools::add.species.to.genus(tree = tree_update_genus, 
                                                      species = insert_genus[i])
  pb_congeneric$tick()
}

# insert in the table
data_insertion[match(insert_genus, data_insertion$s), "insertion"] <- "genus_insertion"

data_round_2 <- data_genus[-which(is.na(match(data_genus$s, tree_update_genus$tip.label)) != TRUE), ]
families_2 <- unique(data_round_2$f)
families_tree <- families_2[!is.na(match(families_2, tree_update_genus$node.label))]
families_not_tree <- families_2[is.na(match(families_2, tree_update_genus$node.label))]
data_round_2 <- data_round_2[data_round_2$f %in% families_tree, ] # only species with data in tree

# family insertion --------------------------------------------------------

tree_update_family <- tree_update_genus
is.ultrametric(tree_update_family) # check here if tree is ultrametric
# phytools::force.ultrametric(tree = tree_update_family, method = "extend")

node_number_family <- lapply(data_round_2$f, function(x) which(c(tree_update_family$tip.label, 
                                          tree_update_family$node.label) == x)) # finding the names of families in tree

pb_family <- progress::progress_bar$new(format = "Adding species to family nodes [:bar] :percent", 
                                            total = length(node_number_family), clear = FALSE, 
                                            width = 60, current = "<", incomplete = ">", 
                                            complete = ">")

node_number_family <- unlist(node_number_family)

for(i in 1:length(node_number_family)){
  tree_update_family <- phytools::bind.tip(tree = tree_update_family, 
                                    tip.label = data_round_2$s[i], where = node_number_family[i], 
                                    position = 0) # position can be a value provided by the user
  # identify if there is any species from the same genus
  pb_family$tick()
}

# species inserted in family level
data_insertion[match(data_round_2$s, data_insertion$s), "insertion"] <- "family_insertion"


## using while to insert since the vector will change the size in non-regular way at each step
#
#spp_to_add_round2 <- data_round_2$s
#pb_family <- progress::progress_bar$new(format = "Adding species to family nodes [:bar] :percent", 
#                                        total = length(spp_to_add_round2), clear = FALSE, 
#                                        width = 60, current = "<", incomplete = ">", 
#                                        complete = ">")
#
#count <- 0
#species_to_genus2 <- vector(mode = "list")
#while (length(spp_to_add_round2) >= 1) {
#  count <- count + 1
#  
#  # family name in which species will be grafted
#  family_name <- data_round_2[match(spp_to_add_round2[1], 
#                              data_round_2$s), "f"]
#  # node family in which species will be grafted
#  node_family <- which(c(tree_update_family$tip.label, 
#                         tree_update_family$node.label) == family_name)
#  
#  tree_update_family <- phytools::bind.tip(tree = tree_update_family, 
#                                           tip.label = spp_to_add_round2[1],
#                                           where = node_family, 
#                                           position = 0)
#  spp_data <- 1:length(spp_to_add_round2)
#  names(spp_data) <- spp_to_add_round2
#  insert_spp <- treedata_modif(phy = tree_update_family, 
#                               data = spp_data, warnings = F)$nc$data_not_tree # refreshing species to be grafted
#  
#  genus_data <- sub("_.*", "", insert_spp)
#  genus_tree <- sub("_.*", "", tree_update_family$tip.label)
#  genus_in_tree <- data_round_2[genus_data %in% genus_tree, "s"]
#  
#  if(is.ultrametric(tree_update_family) != TRUE){
#    tree_update_family <- phytools::force.ultrametric(tree = tree_update_family, method = "extend")
#  }
#  
#  if (length(genus_in_tree) >= 1) {
#    species_to_genus <- genus_in_tree
#    species_to_genus2[[count]] <- species_to_genus
#    for (i in 1:length(species_to_genus)) {
#      tree_update_family <- phytools::add.species.to.genus(tree = tree_update_family, 
#                                                    species = species_to_genus[i])
#    }
#    insert_spp <- treedata_modif(phy = tree_update_family, 
#                                 data = spp_data, warnings = F)$nc$data_not_tree
#  }
#  spp_to_add_round2 <- insert_spp
#  if(is.ultrametric(tree_update_family) != TRUE){
#    tree_update_family <- phytools::force.ultrametric(tree = tree_update_family, method = "extend")
#  }
#  
#  pb_family$tick()
#  
#}

# insertion orders --------------------------------------------------------

tree_update_order <- tree_update_family

data_round_3 <- data[-which(data$s %in% tree_update_order$tip.label == TRUE), ]

# order species
orders <-  tax[which(tax$rank == "order"), ]
orders_no_incertae <- orders[-grep("Incertae sedis in", orders$name), ] # removing incertae sedis
spp_order <- lapply(orders_no_incertae$name, function(x) fishtree::fishtree_taxonomy(ranks = x)[[1]]$sampled_species)
names(spp_order) <- orders_no_incertae$name
spp_order <- lapply(spp_order, function(x) gsub(" ", "_", x))


# naming orders 

not_samp_orders <- names(unlist(lapply(spp_order, function(x) which(length(x) == 0)))) # fdmilies not sampled in phylogeny
monotip_orders <-  names(unlist(lapply(spp_order, function(x) which(length(x) == 1)))) # monotipic families in the tree
one_spp <- c(not_samp_orders, monotip_orders)
order_names_tree <- spp_order[-match(one_spp, names(spp_order))]
position_order <- 
  lapply(order_names_tree, function(x){
    node_anc <- ape::getMRCA(phy = tree_update_order, tip = x)
  })
for(i in 1:length(position_order)){ # naming tree nodes
  tree_update_order$node.label[position_order[[i]] - ape::Ntip(tree_update_order)] <- names(position_order)[i]
}

monotipic_list_order <- spp_order[match(monotip_orders, names(spp_order))]

for(i in 1:length(monotip_orders)){ 
  spp_monotipic <- spp_order[match(monotip_orders[i], names(spp_order))][[1]]
  
  spp_singleton <- tree_update_order$tip.label[match(spp_monotipic, tree_update_order$tip.label)]
  spp_singleton_add <- paste(sub("_.*", "", spp_singleton), 
                             "_", "singleton", sep = "")
  tree_update_order <- phytools::add.species.to.genus(tree = tree_update_order, 
                                                species = spp_singleton_add)
  monotipic_list_order[[i]] <- c(spp_singleton, spp_singleton_add)
}

# naming monotipic orders

position_monotipic <- 
  lapply(monotipic_list_order, function(x){
    node_anc <- ape::getMRCA(phy = tree_update_order, tip = x)
  })
for(i in 1:length(position_monotipic)){ # naming tree nodes
  tree_update_order$node.label[position_monotipic[[i]] - ape::Ntip(tree_update_order)] <- 
    names(position_monotipic)[i]
}

# adding species in order node

node_number_order <- lapply(data_round_3$o, function(x) which(c(tree_update_order$tip.label, 
                                                                 tree_update_order$node.label) == x)) # finding the names of orders in tree
pos_no_order <- which(lapply(node_number_order, function(x) length(x) == 0) == TRUE) # species without order in tree

data_round_3_order <- data_round_3[-pos_no_order, ] # removing species with no order in the tree
node_number_order_tree <- node_number_order[-pos_no_order] # removing orders with no species
pb_order <- progress::progress_bar$new(format = "Adding species to order nodes [:bar] :percent", 
                                       total = length(node_number_order_tree), clear = FALSE, 
                                       width = 60, current = "<", incomplete = ">", 
                                       complete = ">")

for(i in 1:length(node_number_order_tree)){
  tree_update_order <-
    phytools::bind.tip(tree = tree_update_order, 
                       tip.label = data_round_3_order$s[i], 
                       where = node_number_order_tree[i], 
                       position = 0) # position can be a value provided by the user
  # identify if there is any species from the same genus
  pb_order$tick()
}


# final tree and insertion table ------------------------------------------

# croping the tree
data_final_inserted <- data$s[which(data$s %in% tree_update_order$tip.label == TRUE)]
tree_final <- keep.tip(phy = tree_update_order, tip = data_final_inserted)

# insertion table

data_insertion[match(data_round_3_order$s, data_insertion$s), "insertion"] <- "order_insertion"
data_insertion[which(is.na(data_insertion$insertion) == TRUE), "insertion"] <- "not_inserted"
# saving the trees and table --------------------------------------------------------

saveRDS(tree_update_order, here::here("tree_graft_final.rds"))
saveRDS(data_insertion, here::here("data_insertion.rds"))

(table(data_insertion$insertion))/(sum(table(data_insertion$insertion))) * 100
