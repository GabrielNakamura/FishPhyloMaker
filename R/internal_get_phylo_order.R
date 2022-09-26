
get_phylo_order <- function(data){
  rank_order <- as.character(unique(data$o))
  rank_family <- as.character(unique(data$f))
  spp <- as.character(data$s)
  all_families <- unique(unlist(lapply(rank_order, function(x) {
    fishbasedata[which(x == fishbasedata$Order), 5]
  })))
  families_in_orders <- suppressWarnings(all_families[which(unique(data$f) != 
                                                              all_families)])
  families_order_and_data <- unique(c(rank_family, families_in_orders))
  list_family <- vector(mode = "list", length = length(families_order_and_data))
  list_family <- fishtree::fishtree_taxonomy(ranks = families_order_and_data)
  list_family <- lapply(list_family, function(x) x$sampled_species)
  list_family <- lapply(list_family, function(x) gsub(" ", "_", x))
  
  
  monotipic_family <- names(unlist(lapply(list_family, 
                                          function(x) which(length(x) == 1))
  )
  )
  list_monotipic <- vector(mode = "list", length = length(monotipic_family))
  for (i in 1:length(monotipic_family)) {
    list_monotipic[[i]] <- tryCatch(fishtree::fishtree_taxonomy(rank = monotipic_family[i])[[1]]$taxonomy[[9]], 
                                    error = function(e) paste("not.found", "_", monotipic_family[i], 
                                                              sep = ""))
  }
  orders_to_add <- unique(unlist(list_monotipic[-which(sub("_.*", 
                                                           "", unlist(list_monotipic)) == "not.found")]))
  differences_orders_toadd <- setdiff(rank_order, orders_to_add)
  if (length(differences_orders_toadd) >= 1) {
    all_orders_include <- c(differences_orders_toadd, 
                            orders_to_add)
  }
  all_orders_include <- unique(c(rank_order, unique(orders_to_add)))
  list_order <- fishtree::fishtree_taxonomy(ranks = all_orders_include)
  list_order <- lapply(list_order, function(x) x$sampled_species)
  list_order <- lapply(list_order, function(x) gsub(" ", "_", x))
  list_res <- vector(mode = "list", length = 2)
  list_res$family <- list_family
  list_res$order <- list_order
  return(list_res)
}

