# naming node according to order
for (i in 1:length(list_order)) {
  # i = 5
  phylo_family<- ape::makeNodeLabel(phylo_family, "u", nodeList = list(Ord_name = list_order[[i]]))
  phylo_family$node.label[which(phylo_family$node.label == "Ord_name")] <- names(list_order)[i]
}

# naming node according to families 
for (i in 1:length(list_family)) {
  #i= 225
  phylo_family <- ape::makeNodeLabel(phylo_family, "u", nodeList = list(Fam_name = list_family[[i]]))
  phylo_family$node.label[which(phylo_family$node.label == "Fam_name")] <- paste(names(list_family)[i])
}

phylo_family[!is.na(match(list_family[[i]], phylo_family$tip.label))]

families_in_tree <- families_order_and_data[which(!is.na(match(families_order_and_data, phylo_family$node.label)) == T)]
families_monotipic_notfound <- setdiff(monotipic_family, families_in_tree)

for(i in 1:length(families_monotipic_notfound)){
  # i = 73
  
  spp_tmp <-  tryCatch(fishtree::fishtree_taxonomy(rank = families_monotipic_notfound[i])[[1]]$species,
                       error = function(e) paste("not.found", "_", families_monotipic_notfound[i], sep = ""))
  spp_tmp <- gsub("\\ ", "_", spp_tmp)
  
  list_family[families_monotipic_notfound[i]] <- list(spp_tmp)
  
  # list_monotipic[[i]]<- tryCatch(fishtree::fishtree_taxonomy(rank = families_monotipic_notfound[i])[[1]]$species,
  #                                error = function(e) paste("not.found", "_", monotipic_family[i], sep = ""))
}