#' Generate a list of species 
#' Auxiliar function to obtain taxonomic classification and check the names of species
#'   present in species pool
#' 
#' @param data A character vector with species names or a community matrix with species names in columns
#'
#' @return Data frame with three columns containing the name of species (s), the Family (f) and Order (o) 
#' 
#' @export
#'
#' @examples 
#'  \dontrun{
#'     data(neotropical_comm)
#'     data_comm <- neotropical_comm[, -c(1, 2)]
#'     taxon_data <- FishTaxaMaker(data_comm) # formating data
#'     Loricariidae # informing Family of Curculionichthys insperatus
#'     Siluriformes # informing Order of Curculionichthys insperatus
#' }
#' 
#' 
FishTaxaMaker <- function (data) 
{
  if (is.data.frame(data) == TRUE) {
    names_data <- colnames(data)
  }
  if (is.matrix(data) == TRUE) {
    if (dim(data)[2] < 1) {
      stop("\n More than one species must be supplied in occurence matrix \n")
    }
    names_data <- colnames(data)
  }
  if (is.vector(data) == TRUE) {
    if (length(data) == 1) {
      stop("\n More than one species must be supplied \n")
    }
    names_data <- data
  }
  requireNamespace("rfishbase")
  utils::data(fishbase, package = "rfishbase")
  list_genus <- fishbase[match(sub("_.*", "", names_data), 
                               fishbase$Genus), c("Genus", "Family", "Order")]
  list_local <- data.frame(s = names_data, f = list_genus$Family, 
                           o = list_genus$Order)
  not_found <- list_local[which(rowSums(is.na(list_local)) > 
                                  0), ]
  print_cat_Family <- function(not_found_fishtree) {
    cat("tell the Family of ", not_found_fishtree)
    cat("\n")
  }
  print_cat_Order <- function(not_found_fishtree) {
    cat("tell the Order of ", not_found_fishtree)
    cat("\n")
  }
  if (dim(not_found)[1] == 0) {
    sub_Perciformes <- which(list_local$f == "Cichlidae")
    if(length(sub_Perciformes) >= 1){
      list_local[sub_Perciformes, "o"] <- "Cichliformes"
    }
    return(list_local)
  }
  else {
    for (i in 1:length(not_found$s)) {
      spp_family <- readline(prompt = print_cat_Family(not_found_fishtree = not_found$s[i]))
      spp_order <- readline(prompt = print_cat_Order(not_found_fishtree = not_found$s[i]))
      list_local[match(not_found$s[i], list_local$s), 
                 "f"] <- spp_family
      list_local[match(not_found$s[i], list_local$s), 
                 "o"] <- spp_order
    }
    sub_Perciformes <- which(list_local$f == "Cichlidae")
    if(length(sub_Perciformes) >= 1){
      list_local[sub_Perciformes, "o"] <- "Cichliformes"
    }
    return(list_local)
  }
  if(dim(not_found)[1] >= 1){
    warning(paste("Species", not_found, "was not found in Fishbase", sep = ""))
  }
}