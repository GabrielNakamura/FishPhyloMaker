#' Internal function to help in interactive process 
#'
#' @param print_cat Character
#' @param spp Character
#' @param order Character
#'
#' @return interactive action in consele
#'
print_cat_family <- function(print_cat, spp, order){
  cat("in which/between Family species", spp, "must be inserted")
  cat("\n")
  for(i in 1:length(print_cat)){
    cat(i, print_cat[i], sep="\t")
    cat("\n")
  }
  cat("\n")
  cat(i + 1, "or, to insert in the node corresponding to the Order type:", order, sep= "\t")
  cat("\n")
}