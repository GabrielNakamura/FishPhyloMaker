genus_insertion <- function(insert_spp, phy){
  
  if (progress.bar == TRUE) {
    pb_congeneric <- progress::progress_bar$new(format = "Adding congeneric species [:bar] :percent", 
                                                total = length(insert_spp), clear = FALSE, 
                                                width = 60, current = "<", incomplete = ">", 
                                                complete = ">")
  }
  
  for (i in 1:length(insert_spp)) {
    phy <- phytools::add.species.to.genus(tree = phy, 
                                                  species = insert_spp[i])
    if (progress.bar == TRUE) {
      pb_congeneric$tick()
    }
  }
  return(phy)
}