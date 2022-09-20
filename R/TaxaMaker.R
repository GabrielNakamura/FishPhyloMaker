test_fish_list = tibble::tibble(
  species = c("Serrasalmus_geryi", "Careproctus_reinhardti", "Gobiomorphus_coxii", 
              "Periophthalmus_barbarus", "Prognichthys_glaphyrae", "Barathronus_bicolor", 
              "Knipowitschia_croatica", "Rhamphochromis_lucius", "Neolissochilus_tweediei", 
              "Haplochromis_nyanzae", "Astronesthes_micropogon", "Sanopus_reticulatus"),
  genus = c("Serrasalmus", "Careproctus", "Gobiomorphus", "Periophthalmus",
            "Prognichthys", "Barathronus", "Knipowitschia", "Rhamphochromis", 
            "Neolissochilus", "Haplochromis", "Astronesthes", "Sanopus"),
  family = c("Serrasalmidae", "Liparidae", "Eleotridae", "Gobiidae", 
             "Exocoetidae", "Aphyonidae", "Gobiidae", "Cichlidae", 
             "Cyprinidae", "Cichlidae", "Stomiidae", "Batrachoididae")
)

data = test_fish_list$species
allow.manual.insert = TRUE
user_taxon_data = NULL

FishTaxaMaker <- function (data, 
                           allow.manual.insert = TRUE) # look for checking databases for taxonomic names 
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
  
  l_cas <- vector(mode = "list", length = length(names_data)) 
  pb = txtProgressBar(min = 0, max = length(l_cas), initial = 1, style = 3) 
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
  pos <- which(vapply(l_cas, function(x) is.null(x), TRUE) == TRUE) 
  l_cas[pos] <- paste("not_found", names_data[pos], sep = "_") 
  names_all <- 
  rep.int(names_data, times = unlist(lapply(l_cas, function(x){
    dim(data.frame(x))[1]
  } )))
  df_cas <- do.call(rbind, l_cas)
  df_cas$name.user <- names_all # queries from catalog of fishes

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
  # joining both data
  df_join[grep("not_found", df_join[, "species"]), "species"] <- NA
  all_fb_names <- gsub(" ", "_", c(df_taxon$valid_names, not_found_fishtree))
  fb_names_org <- table(df_cas$name.user)[match(all_fb_names, names(table(df_cas$name.user)))]
  names_fb_join <- rep.int(names(fb_names_org), times = fb_names_org)
  df_join <- data.frame(df_cas, fishbase = names_fb_join) # both fishbase and cas
  df_join <- df_join[, c("name.user", "species", "fishbase", "status", "family", "query")]
  df_join$species <- gsub(" ", "_", df_join$species)
  
  consolidate <- apply(df_join, 1, function(x) {
    df_join$consolidate <- NA
    df_join[df_join$species %in% df_join$fishbase, "consolidate"] <- df_join[df_join$species %in% df_join$fishbase, "species"]
    df_join[which(is.na(df_join$species) & is.na(df_join$fishbase) == TRUE), "consolidate"] <- NA # no valid names
    df_join[which(is.na(df_join$species) & !is.na(df_join$fishbase) == TRUE), "consolidate"] <- df_join[which(is.na(df_join$species) & !is.na(df_join$fishbase) == TRUE), "fishbase"] # species only fb
    df_join[which(!is.na(df_join$species) & is.na(df_join$fishbase) == TRUE), "consolidate"] <- df_join[which(!is.na(df_join$species) & is.na(df_join$fishbase) == TRUE), "species"] # species only CAS
    df_join[which(!is.na(df_join$species) & !is.na(df_join$fishbase) == TRUE), "consolidate"] <- df_join[which(!is.na(df_join$species) & !is.na(df_join$fishbase) == TRUE), "species"] # species in both
  })
  
  
  list_res <- vector(mode = "list", length = 3)
  tax_hierarch <- fishbasedata[match(df_taxon$valid_names, 
                                     fishbasedata$Species), c("Subfamily", "Family", "Order", 
                                                              "Class", "SuperClass")]
  data_fishbase_complete <- cbind(df_taxon, tax_hierarch)
  data_fishbase_complete$valid_names <- gsub(" ", "_", data_fishbase_complete$valid_names)
  list_res[[1]] <- data_fishbase_complete
  list_res[[2]] <- data_fishbase_complete[, c("valid_names", 
                                              "Family", "Order")]
  colnames(list_res[[2]]) <- c("s", "f", "o")
  list_res[[2]] <- list_res[[2]][match(unique(list_res[[2]]$s), 
                                       list_res[[2]]$s), ]
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
