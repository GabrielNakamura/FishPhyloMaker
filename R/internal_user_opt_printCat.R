print_cat<- function(print_cat, spp){
  cat("in which species species", spp, "must be inserted")
  cat("\n")
  for(i in 1:length(print_cat)){
    cat(i, print_cat[i], sep="\t")
    cat("\n")
  }
  cat("\n")
  cat(i+1, "family_node", sep= "\t")
  cat("\n")
}