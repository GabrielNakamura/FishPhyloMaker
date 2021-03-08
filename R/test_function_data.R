library(fishtree)
library(rfishbase)
library(phytools)

data_all <- load(here::here("data", "neotropical_comm.rda"))
data_comm <- neotropical_comm[, -c(1, 2)]

source(here::here("R", "tab_function.R"))
source(here::here("R", "internal_user_opt_printCat.R"))
source(here::here("R", "internal_user_opt_printCat2.R"))
source(here::here("R", "internal_user_opt_printCatFamily.R"))
source(here::here("R", "internal_treedata_modif.R"))
source(here::here("R", "internal_filter_rank.R"))
source(here::here("R", "test_major_changes.R"))
taxon_data <- tab_function(data_comm)
Loricariidae
Siluriformes

data_process <- taxon_data[c(1, 2, 6, 29, 58), ]
data_process <- rbind(data_process, c("Curculionichthys_inesperado", "Loricariidae", "Siluriformes") )
data_process <- rbind(data_process, c("Dinkiwinki_dipsii", "Teletubiidae", "Cichliformes"))
data_process <- rbind(data_process, c("Peixe_loricariaentregeneros", "Loricariidae", "Siluriformes"))
data_process <- rbind(data_process, c("Peixo_basefamilia", "Loricariidae", "Siluriformes"))
data <- data_process

res_test <- FishPhyloMaker(data = data, return.insertions = TRUE)
quartz()
plot(res_test$Phylogeny, show.node.label = T)

