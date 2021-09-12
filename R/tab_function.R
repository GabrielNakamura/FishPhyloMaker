#' Generate a list of species 
#' Auxiliary function to obtain taxonomic classification and check the names of species
#'   present in species pool
#' 
#' @param data A character vector with species names or a community matrix with species names in columns
#' @param allow.manual.insert Logical, if TRUE (default), the user must type the names of Family and Order of species
#'    not found in Fishbase
#'
#' @return List with three elements. 
#' 
#'     - A data frame containing the taxonomic classification of valid species accordingy to Fishbase
#'     
#'     - A data frame with three columns containing the name of species (s), the Family (f) and Order (o) that can be used in
#'        FishPhyloMaker function
#'        
#'     - A character vector containing all names of species that was not find in Fishbase
#' 
#' @importFrom stats na.omit
#' 
#' @export
#'
#' @examples 
#'  \dontrun{
#'  data(neotropical_comm)
#'  data_comm <- neotropical_comm[, -c(1, 2)]
#'  taxon_data <- FishTaxaMaker(data_comm, allow.manual.insert = TRUE)
#'  Characidae
#'  Characiformes
#'  Characidae
#'  Characiformes
#'  Characidae
#'  Characiformes
#'  Loricariidae
#'  Siluriformes
#'  Characidae
#'  Characiformes
#'  Cichlidae
#'  Cichliformes
#'  Crenuchidae
#'  Characiformes
#'  Gymnotidae
#'  Gymnotiformes
#'  Loricariidae
#'  Siluriformes
#'  Loricariidae
#'  Siluriformes
#'  Loricariidae
#'  Siluriformes
#'  Loricariidae
#'  Siluriformes
#'  Heptapteridae
#'  Siluriformes
#'  Characidae
#'  Characiformes
#'  Loricariidae
#'  Siluriformes
#'  Characidae
#'  Characiformes
#' }
#' 
#'
FishTaxaMaker <- function (data, allow.manual.insert = TRUE) 
{
  if (is.data.frame(data) == TRUE) {
    names_data <- colnames(data)
    names_data_rfishbase <- gsub("_", " ", names_data)
  }
  if (is.matrix(data) == TRUE) {
    if (dim(data)[2] < 1) {
      stop("\n More than one species must be supplied in occurence matrix \n")
    }
    names_data <- colnames(data)
    names_data_rfishbase <- gsub("_", " ", names_data)
  }
  if (is.vector(data) == TRUE) {
    if (length(data) == 1) {
      stop("\n More than one species must be supplied \n")
    }
    names_data <- data
    names_data_rfishbase <- gsub("_", " ", names_data)
  }
  df_taxon <- data.frame(user_spp = names_data, valid_names = rep(NA, 
                                                                  length(names_data)))
  fishbasedata <- as.data.frame(data.frame(rfishbase::load_taxa()))
  not_find <- names_data_rfishbase[which(is.na(match(names_data_rfishbase, 
                                                     fishbasedata$Species)) == TRUE)]
  found <- gsub("_", " ", names_data[which(!is.na(match(names_data_rfishbase, 
                                                        fishbasedata$Species)) == TRUE)])
  names_not_find <- rfishbase::synonyms(species_list = not_find)
  names_not_find_valid <- names_not_find[which(names_not_find$Status == 
                                                 "synonym"), ]
  df_taxon[match(found, gsub("_", " ", df_taxon$user_spp)), 
           "valid_names"] <- found
  df_taxon[match(names_not_find_valid$synonym, gsub("_", " ", 
                                                    df_taxon$user_spp)), "valid_names"] <- names_not_find_valid$Species
  not_found_fishtree <- data.frame(names_not_find[match(gsub("_", 
                                                             " ", df_taxon[which(is.na(df_taxon$valid_names) == TRUE), 
                                                                           "user_spp"]), names_not_find$synonym), ])$synonym
  list_res <- vector(mode = "list", length = 3)
  tax_hierarch <- fishbasedata[match(df_taxon$valid_names, 
                                     fishbasedata$Species), c("Subfamily", "Family", "Order", 
                                                              "Class", "SuperClass")]
  data_fishbase_complete <- cbind(df_taxon, tax_hierarch)
  list_res[[1]] <- data_fishbase_complete
  list_res[[2]] <- data_fishbase_complete[, c("valid_names", 
                                              "Family", "Order")]
  colnames(list_res[[2]]) <- c("s", "f", "o")
  list_res[[2]] <- list_res[[2]][match(unique(list_res[[2]]$s), 
                                       list_res[[2]]$s), ]
  list_res[[2]]$s <- gsub(" ", "_", list_res[[2]]$s)
  list_res[[2]][match(gsub(" ", "_", na.omit(not_found_fishtree)), 
                      list_res[[1]]$user_spp), "s"] <- na.omit(not_found_fishtree)
  list_res[[2]][match(gsub(" ", "_", na.omit(not_found_fishtree)), 
                      list_res[[1]]$user_spp), c("f", "o")] <- paste("not_find")
  list_res[[2]]$s <- gsub(" ", "_", list_res[[2]]$s)
  if (length(not_found_fishtree) >= 1) {
    list_res[[3]] <- not_found_fishtree
  }
  else {
    list_res[[3]] <- paste("All species were found in Fishtree")
  }
  names(list_res) <- c("All_info_fishbase", "Taxon_data_FishPhyloMaker", 
                       "Species_not_in_Fishbase")
  if (length(not_found_fishtree) >= 1) {
    if (allow.manual.insert == TRUE) {
      print_cat_Family <- function(not_found_fishtree) {
        cat("tell the Family of ", not_found_fishtree)
        cat("\n")
      }
      print_cat_Order <- function(not_found_fishtree) {
        cat("tell the Order of ", not_found_fishtree)
        cat("\n")
      }
      for (i in 1:length(not_found_fishtree)) {
        spp_family <- readline(prompt = print_cat_Family(not_found_fishtree = not_found_fishtree[i]))
        spp_order <- readline(prompt = print_cat_Order(not_found_fishtree = not_found_fishtree[i]))
        list_res[[2]][which(list_res[[1]]$user_spp == 
                              gsub(" ", "_", not_found_fishtree[i])), c("s", 
                                                                        "f")] <- c(not_found_fishtree[i], spp_family)
        list_res[[2]][which(list_res[[2]]$s == not_found_fishtree[i]), 
                      "o"] <- spp_order
      }
      list_res[[2]]$s <- gsub(" ", "_", list_res[[2]]$s)
    }
  }
  return(list_res)
}
