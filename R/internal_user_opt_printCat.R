#' Internal function - auxiliar to interactive procedure
#'
#' @param print_cat Character
#' @param spp Character
#' @param family Character
#'
#' @return interactive action in console
#' @keywords internal
#' 
print_cat <- function(print_cat, spp, family){
  cat("in which/between Genus species", spp, "must be inserted")
  cat("\n")
  for(i in 1:length(print_cat)){
    cat(i, print_cat[i], sep="\t")
    cat("\n")
  }
  cat("\n")
  cat(i + 1, "or, to insert in the node corresponding to the Family:", family, sep= "\t")
  cat("\n")
}