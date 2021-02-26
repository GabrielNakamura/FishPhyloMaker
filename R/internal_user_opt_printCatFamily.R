#' Internal function to help in interactive process 
#'
#' @param print_cat Character
#' @param spp Character
#' @param family Character
#'
#' @return
#'
#' @examples
print_cat_family <- function(print_cat, spp, family){
  cat("in which/between Family species", spp, "must be inserted")
  cat("\n")
  for(i in 1:length(print_cat)){
    cat(i, print_cat[i], sep="\t")
    cat("\n")
  }
  cat("\n")
  cat(i + 1, "or, to insert in the node corresponding to the Order type:", family, sep= "\t")
  cat("\n")
}