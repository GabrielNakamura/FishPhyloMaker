
get_phylo_order <- function(data){
  rank_order <- as.character(unique(data$o))
  rank_family <- as.character(unique(data$f))
  spp <- as.character(data$s)
  all_families <- unique(fishbasedata[fishbasedata$Order %in% data$o, "Family"]) # all families in orders in data
  families_order_and_data <- unique(c(data$f, all_families))
  list_family <- fishtree::fishtree_taxonomy(ranks = families_order_and_data) 
  list_family <- lapply(list_family, function(x) x$sampled_species)
  list_family <- lapply(list_family, function(x) gsub(" ", "_", x))
  not_found_families <- c(setdiff(families_order_and_data, names(list_family)))
  monotipic_family <- names(unlist(lapply(list_family, function(x) which(length(x) <= 1)))) # monotipic families in tree
  monotipic_family <- c(not_found_families, monotipic_family) # monotipic and families not sampled in tree
  list_monotipic <- 
  lapply(monotipic_family, function(x){
    tryCatch(fishtree::fishtree_taxonomy(rank = x)[[1]]$taxonomy[[9]], 
             error = function(e) paste("not.found", "_", x, 
                                       sep = ""))
  })
  
  orders_to_add <- unique(unlist(list_monotipic[-which(sub("_.*", 
                                                           "", unlist(list_monotipic)) == "not.found")])) # orders with sampled families in tree
  differences_orders_toadd <- setdiff(rank_order, orders_to_add)
  if (length(differences_orders_toadd) >= 1) { # all orders 
    all_orders_include <- c(differences_orders_toadd, 
                            orders_to_add)
  }
  all_orders_include <- unique(c(rank_order, unique(orders_to_add)))
  list_order <- lapply(all_orders_include, function(x){
    tryCatch(fishtree::fishtree_taxonomy(ranks = x)[[1]]$sampled_species, 
             error = function(e) paste("not.found", "_", x, 
                                       sep = ""))
  })
  names(list_order) <- all_orders_include
  list_order <- lapply(list_order, function(x) gsub(" ", "_", x))
  list_res <- vector(mode = "list", length = 5)
  list_res$family <- list_family
  list_res$order <- list_order
  list_res$families_order_and_data <- families_order_and_data
  list_res$not_found_families <- not_found_families
  list_res$monotipic_families <- monotipic_family
  return(list_res)
}
