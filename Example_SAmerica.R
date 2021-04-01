devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main")

library(FishPhyloMaker)
data("fish_SAmerica")

taxon_data <- FishTaxaMaker(data = fish_SAmerica)
Loricariidae
Siluriformes
Anostomidae 
Characiformes
Anostomidae 
Characiformes
Anostomidae 
Characiformes
Serrasalmidae
Characiformes
Auchenipteridae
Siluriformes
Auchenipteridae
Siluriformes
Auchenipteridae
Siluriformes
str(taxon_data)

phylo_fish_SAmerica <- FishPhyloMaker(data = taxon_data, return.insertions = T)
Serrasalmidae
Acanthicus
Schizodon
Liosomadoras
Hemiancistrus
Ctenolucius
Roeboexodon Exodon
Leporinus
Gymnotus
Loricariichthys
Anostomus
Furcodontichthys
Lasiancistrus
Leporinus
Schizodon Hypomasticus
Batrochoglanis
Pseudolithoxus
Leporinus
Anodus
Paralonchurus
Pachypops
Vandellia Paravandellia
Pachyurus
Anostomus
Sternarchorhamphus
Vandellia
Odontognathus
Exallodontus
Myleus
Gnathodolus
Trachelyopterus
Hemiodontichthys
Rhinelepis
Pseudohemiodon
Schizodon
Gnathodolus
Panaque Panaqolus
Pseudacanthicus Leporacanthicus
Hypopomidae
Sternarchogiton
Platyurosternarchus
Gymnocorymbus
Schizodon
Synbranchidae
Ageneiosus
Asterophysus
Trachelyopterus
Ageneiosus

library(phytools)
plotTree(phylo_fish_SAmerica$Phylogeny, fsize = 0.1, ftype = "i",
         type = "fan", lwd = 0.3)

tree <- phylo_fish_SAmerica$Phylogeny
tree <- ape::makeNodeLabel(tree)

phylo <- tree

rm.famNames <- which(table(taxon_data$f) == 1) # monotipic families
names.fam <- setdiff(unique(taxon_data$f), names(rm.famNames))

for (i in 1:length(names.fam)) {
  #i=3
  set <- subset(taxon_data, f == names.fam[i])
  phylo <- ape::makeNodeLabel(phylo, "u", nodeList = list(Fam_name = set$s))
  
  phylo$node.label[which(phylo$node.label == 
                           "Fam_name") ] <- paste(set$f[1])
}


pos.node <- unlist(lapply(names.fam, function(x){
  which(phylo$node.label == x) + length(tree$tip.label)
}))

df.phylo <- data.frame(Fam.names = names.fam[-33],
                       node.number = pos.node[-33])

library(viridis)
p.base <- ggtree(phylo, layout = "circular") + # geom_tiplab2(size = 0.4) +
  geom_hilight(data = df.phylo, aes(node = node.number, fill = Fam.names), 
               alpha = .6) +
  scale_fill_viridis(discrete = T) + theme(legend.position = "none")
  
p.base
