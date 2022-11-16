library(FishPhyloMaker)
t <- readRDS("/Users/gabrielnakamura/OneDrive/Manuscritos/Darwinian_shortfalls/output/data_insertion.rds")
tree <- readRDS("/Users/gabrielnakamura/OneDrive/Manuscritos/Darwinian_shortfalls/output/tree_graft_final.rds")
orders <- names(sort(table(t$o), decreasing = TRUE))[1:10] # top ten orders with most species
t_order <- lapply(orders, function(x) t[which(t$o == x), ])  # all specie sfor each one of the ten orders
sub_samp <- seq(3, 100, 10) # sub samples for the bootstrap procedure

lapply(sub_samp, function(x){
  sub_insert <- t_order[[1]][sample(1:nrow(t_order[[1]]), x), ]
  PD_deficit(phylo = tree, data = sub_insert, level = c(""))
})


rep_boot <- rep(sub_samp, each = 10)
boot_1 <- 
lapply(rep_boot, function(x){
  sub_insert <- t_order[[1]][sample(1:nrow(t_order[[1]]), x), ]
  PD_deficit_2(phylo = tree,
               data = sub_insert, 
               level = c("genus_insertion",
                         "family_insertion", 
                         "order_insertion")
  )
})
boot_deficit <- do.call(rbind, boot_1)
names_boot <- paste("boot", rep_boot, sep = "_")
data.frame(boot_deficit, names_boot)
sub_insert <- t_order[[1]][sample(1:nrow(t_order[[1]]), 3), ]
PD_deficit_2(phylo = tree,
           data = sub_insert, 
           level = c("genus_insertion",
                     "family_insertion", 
                     "order_insertion")
           )
PD_deficit_2(phylo = res_phylo_marine$Phylogeny, 
           data = insertions_marine, 
           level = c("genus_insertion", 
                     "family_insertion",
                     "order_insertion")
)
PD_curves_orders <- 
lapply(t_order, function(x){
  PD_defict_curves(phy = tree, data = x, reps = 10, samp = sub_samp)
})
all_deficit <- do.call(rbind, PD_curves_orders[2:10])
all_deficit$group <- rep(paste("order", c(2:10)), each = 1500)
all_deficit$names_2 <- factor(all_deficit$names_boot, levels = unique(all_deficit$names_boot))


ggplot(data = all_deficit, aes(x = names_2, y = Darwinian_deficit, group = group, color = group)) +
  geom_smooth()

ggplot() +
  geom_point(data = mean_deficit, aes(x = names_2, y = mean)) +
  geom_line()

match(PD_curves_orders[[1]]$names_boot)
PD_defict_curves(phy = tree, data = sub_insert, reps = 50, samp = sub_samp)