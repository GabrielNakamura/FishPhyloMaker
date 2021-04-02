devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main", force = T)

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

################################################################

library(ggtree)
library(ggplot2)
library(viridis)

# plotting the phylogenetic tree, accounting with the species insertions types
tree <- phylo_fish_SAmerica$Phylogeny

# finding only the species that was add in the phylogenetic tree
insertions <- phylo_fish_SAmerica$Insertions_data[-which("Present_in_Tree" == 
                                                            phylo_fish_SAmerica$Insertions_data[, "insertions"]), ]

p.base <- ggtree(tree, layout = "circular", size = .3)  %<+% insertions +
  geom_treescale(x = 0, width = 20, linesize = .5, color = "blue", 
                fontsize = 0) + #  plot the scale bar
  annotate("text", x = 4, y = 500, label = "20 myr", size = 1.5) # an attempt for add a scale bar

p.full <- p.base +
  geom_tippoint(aes(color = insertions), 
                size = .5, alpha = .8) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 2)))  +
  scale_color_viridis_d(name = NULL, na.translate = F,
                       labels = c("Congeneric", "Congeneric F", "Family")) 


png(here::here("vignettes", "phylo_SAmerica.png"), width = 10, height = 8, 
    units = "cm", res = 300)
p.full
dev.off()

###################################################################
# script for highlight the family nodes

phylo <- ape::makeNodeLabel(tree)

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

p.base <- ggtree(phylo, layout = "circular") + # geom_tiplab2(size = 0.4) +
  geom_hilight(data = df.phylo, aes(node = node.number, fill = Fam.names), 
               alpha = .6) +
  scale_fill_viridis(discrete = T, option = "B") + theme(legend.position = "none")
  
p.base
