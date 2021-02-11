print_cat<- function(print_cat, spp, family){
  cat("in which/between genus species", spp, "must be inserted")
  cat("\n")
  for(i in 1:length(print_cat)){
    cat(i, print_cat[i], sep="\t")
    cat("\n")
  }
  cat("\n")
  cat(i+1, "or, to insert in family_node type:", family, sep= "\t")
  cat("\n")
}