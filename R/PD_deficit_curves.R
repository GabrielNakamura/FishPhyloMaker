
PD_defict_curves <- 
function(phy, data, reps, samp){
  rep_boot <- rep(samp, each = reps)
  boot_1 <- 
    lapply(rep_boot, function(x){
      sub_insert <- data[sample(1:nrow(data), x, replace = TRUE), ]
      PD_deficit_2(phylo = phy,
                   data = sub_insert, 
                   level = c("genus_insertion",
                             "family_insertion", 
                             "order_insertion")
      )
    })
  boot_deficit <- do.call(rbind, boot_1)
  names_boot <- paste("boot", rep_boot, sep = "_")
  res <- data.frame(boot_deficit, names_boot)
  return(res)
}
