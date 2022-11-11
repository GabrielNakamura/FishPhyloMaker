
FishTaxaMaker <- function (data, 
                           allow.manual.insert = FALSE) # look for checking databases for taxonomic names 
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
                                                                  length(names_data))) # data frame to receive status 
  
  l_cas <- vector(mode = "list", length = length(names_data)) # list to receive cas name search
  pb = txtProgressBar(label = "searching data in CAS - Catalog of fishes", min = 0, max = length(l_cas), initial = 1, style = 3) 
  for(i in 1:length(names_data)){ # searching names in catalog of fishes 
    # i = 1
    tryCatch(
      l_cas[[i]]<- rFishTaxa::search_cas(query = names_data[i], type = "species"),
      error = function(e){
       l_cas[[i]] <-  paste("not_found", names_data[i], sep = "_")
       
      }
    )
    setTxtProgressBar(pb = pb, value = i)
  }
  pos <- which(vapply(l_cas, function(x) is.null(x), TRUE) == TRUE) # which species failed during the search 
  l_cas[pos] <- paste("not_found", names_data[pos], sep = "_") # species that failed (probably scrapping error)
  names_all <- 
  rep.int(names_data, times = unlist(lapply(l_cas, function(x){
    dim(data.frame(x))[1]
  } ))) # repeating queried names to match with the length from cas
  df_cas <- do.call(rbind, l_cas) # cas data
  df_cas$name.user <- names_all # joining names provided by user with names from cas

  fishbasedata <- as.data.frame(data.frame(rfishbase::load_taxa())) # all data from fishbase
  not_find <- names_data_rfishbase[which(is.na(match(names_data_rfishbase, 
                                                     fishbasedata$Species)) == TRUE)] # names not found in fishbase
  found <- gsub("_", " ", names_data[which(!is.na(match(names_data_rfishbase, 
                                                        fishbasedata$Species)) == TRUE)]) # species found in fishbase
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
  # joining both data

  all_fb_names <- gsub(" ", "_", c(df_taxon$valid_names, not_found_fishtree))
  fb_names_org <- table(df_cas$name.user)[match(all_fb_names, names(table(df_cas$name.user)))]
  names_fb_join <- rep.int(names(fb_names_org), times = fb_names_org)
  df_join <- data.frame(df_cas, fishbase = names_fb_join) # both fishbase and cas
  df_join <- df_join[, c("name.user", "species", "fishbase", "status", "family", "query")]
  df_join$species <- gsub(" ", "_", df_join$species)
  df_join[grep("not_found", df_join$species), "species"] <- NA # NA to those species that were not found in cas
  
  # consolidation fb and cas
  df_join$consolidate <- NA
  df_join[df_join$species %in% df_join$fishbase, "consolidate"] <- df_join[df_join$species %in% df_join$fishbase, "species"]
  df_join[which(is.na(df_join$species) & is.na(df_join$fishbase) == TRUE), "consolidate"] <- NA # no valid names
  df_join[which(is.na(df_join$species) & !is.na(df_join$fishbase) == TRUE), "consolidate"] <- df_join[which(is.na(df_join$species) & !is.na(df_join$fishbase) == TRUE), "fishbase"] # species only fb
  df_join[which(!is.na(df_join$species) & is.na(df_join$fishbase) == TRUE), "consolidate"] <- df_join[which(!is.na(df_join$species) & is.na(df_join$fishbase) == TRUE), "species"] # species only CAS
  df_join[which(!is.na(df_join$species) & !is.na(df_join$fishbase) == TRUE), "consolidate"] <- df_join[which(!is.na(df_join$species) & !is.na(df_join$fishbase) == TRUE), "species"] # species in both

  list_res <- vector(mode = "list", length = 3)
  tax_hierarch <- fishbasedata[match(gsub("_", " ", df_join$consolidate), 
                                     fishbasedata$Species), c("Subfamily", "Family", "Order", 
                                                              "Class", "SuperClass")] # whole hierarchy from fishbase
  data_fishbase_cas <- cbind(df_join, tax_hierarch) # joining both databases
  data_fishbase_cas <- data_fishbase_cas[, c("name.user", "species", "fishbase", "status", "consolidate", 
                                             "Subfamily", "Family", "Order", "Class", "SuperClass")]
  colnames(data_fishbase_cas) <- c("user.names", "cas", "fishbase", "status", "consolidate",
                                   "Subfamily", "Family", "Order", "Class", "SuperClass")
  
  list_res[[1]] <- data_fishbase_cas
  list_res[[2]] <- data_fishbase_cas[, c("consolidate", 
                                              "Family", "Order")]
  colnames(list_res[[2]]) <- c("s", "f", "o")
  list_res[[2]] <- list_res[[2]][match(unique(list_res[[2]]$s), 
                                       list_res[[2]]$s), ]
  list_res[[2]][match(gsub(" ", "_", na.omit(not_found_fishtree)), 
                      list_res[[1]]$user_spp), "s"] <- na.omit(not_found_fishtree)
  list_res[[2]][match(gsub(" ", "_", na.omit(not_found_fishtree)), 
                      list_res[[1]]$user_spp), c("f", "o")] <- paste("not_found")
  if (length(not_found_fishtree) >= 1) {
    list_res[[3]] <- not_found_fishtree
  }
  else {
    list_res[[3]] <- paste("All species were found in Fishtree")
  }
  names(list_res) <- c("All_fishbase_cas", "Taxon_data_FishPhyloMaker", 
                       "Species_not_in_Fishbase")
  if (length(not_found_fishtree) >= 1) {
    if (allow.manual.insert == TRUE) { # manual insertion procedure
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
        list_res[[2]][which(list_res[[1]]$user_spp == 
                              gsub(" ", "_", not_found_fishtree[i])), c("s", 
                                                                        "f")] <- c(not_found_fishtree[i], spp_family)
        
        spp_order <- unique(fishbasedata[which(fishbasedata$Family == spp_family), "Order"])
        if(length(spp_order) == 0){
          warning("\n Please, check if the family name typed is valid \n")
          spp_order <- readline(prompt = print_cat_Order(not_found_fishtree = not_found_fishtree[i]))
        }
        list_res[[2]][which(list_res[[2]]$s == not_found_fishtree[i]), 
                      "o"] <- spp_order
      }
      list_res[[2]]$s <- gsub(" ", "_", list_res[[2]]$s)
    }
  }
  return(list_res)
}
