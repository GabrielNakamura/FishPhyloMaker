#' Internal function
#'
#' @param spp Character
#'
#' @return interactive action in console
#'
print_cat2<- function(spp){
  cat("To insert the following species:", spp, "provide a newick file or type politomy to insert species as politomy
      in its respective order node")
  cat("\n")
}
