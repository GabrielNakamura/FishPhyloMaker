

# order for monotipic families
monotipic_family <- names(unlist(lapply(list_family, function(x) which(length(x) == 1))))
list_monotipic <- vector(mode = "list", length = length(monotipic_family))
for(i in 1:length(monotipic_family)){
  list_monotipic[[i]]<- tryCatch(fishtree::fishtree_taxonomy(rank = monotipic_family[i])[[1]]$taxonomy[[9]],
                              error = function(e) paste("not.found", "_", monotipic_family[i], sep = ""))
}
orders_to_add <- unique(unlist(list_monotipic[- which(sub("_.*", "", unlist(list_monotipic)) == "not.found")]))
setdiff(rank_order, orders_to_add)
fishtree::fishtree_taxonomy(rank = monotipic_family[22])
list_monotipic[which(sub("_.*", "", unlist(list_monotipic)) == "not.found")]
match("Characidae", phylo_order$node.label)


fishtree::fishtree_phylogeny(rank = "Incertae sedis in Eupercaria")
fishtree::fishtree_taxonomy(rank = "Teletubiidae")[[1]]$taxonomy[[9]]
fishtree::fishtree_taxonomy(rank = "Zanclidae")$Zanclidae$taxonomy$order


which(fishtree_phylogeny(rank = fishtree::fishtree_taxonomy(rank = "Zanclidae")$Zanclidae$taxonomy$order)$tip.label == "Zanclus_cornutus")


match(list_family$Synbranchidae, phylo_order$tip.label)
ape::makeNodeLabel(phylo_order, "u", nodeList = list(Fam_name = list_family$Synbranchidae))

unlist(list_family$Zanclidae[!is.na(match(list_family$Zanclidae, phylo_order$tip.label))])

families_not_found_fishtree <- names(unlist(lapply(lapply(list_family, function(x) {
  sub("_.*", "", x)
}), function(y) which(length(y) == 1) & which(y == "not.found"))))

list_family_tobeaddnames <- list_family[- match(families_not_found_fishtree, names(list_family))]

fishtree::fishtree_phylogeny(rank = "Latidae")
fishtree::fishtree_phylogeny(rank = "Zanclidae")
fishtree::fishtree_phylogeny(rank = "Teletubidae")
fishtree::fishtree_phylogeny(rank = "Scaridae")
fishtree::fishtree_phylogeny(rank = "Acanthuriformes")


fishtree::fishtree_taxonomy(rank = "Pholidichthyiformes")
