data(neotropical_comm)
data_comm <- neotropical_comm[, -c(1, 2)] # removing latitude and longitude

library(FishPhyloMaker)
taxon_dataPR <- FishTaxaMaker(data_comm)
Loricariidae
Siluriformes

res_phylo <- FishPhyloMaker(data = taxon_dataPR, return.insertions = TRUE)
Hisonotus 
Synbranchidae 

tree.PR<- res_phylo$Phylogeny
tree.PR <- ape::makeNodeLabel(tree.PR)
phylo <- tree.PR

rm.famNames <- which(table(taxon_dataPR$f) == 1) # monotipic families
names.fam <- setdiff(unique(taxon_dataPR$f), names(rm.famNames))

for (i in 1:length(names.fam)) {
  #i=3
  set <- subset(taxon_dataPR, f == names.fam[i])
  phylo <- ape::makeNodeLabel(phylo, "u", nodeList = list(Fam_name = set$s))
  
  phylo$node.label[which(phylo$node.label == 
                           "Fam_name") ] <- paste(set$f[1])
}

pos.node <- unlist(lapply(names.fam, function(x){
  which(phylo$node.label == x) + length(phylo$tip.label)
}))

df.phylo <- data.frame(Fam.names = names.fam,
                       node.number = pos.node)

ggtree(phylo, layout = "circular") + geom_tiplab2(size = 2) +
  geom_hilight(data = df.phylo, aes(node = node.number, fill = Fam.names), 
               alpha = .6) +
  scale_fill_viridis(discrete = T) #+ theme(legend.position = "none")


