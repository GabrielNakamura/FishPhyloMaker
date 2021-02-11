#' Generate a list of species 
#'
#' @param data A character vector with species names or a community matrix with species names in columns
#'
#' @return
#' @export
#'
#' @examples
#' 
tab_function<- function(data){
  if(is.data.frame(data) == TRUE){
    names_data <- colnames(data)
  }
  if(is.matrix(data) == TRUE){
    if(dim(data)[2] < 1){
      stop("\n More than one species must be supplied in occurence matrix \n")
    }
    names_data <- colnames(data)
  }
  if(is.vector(data) == TRUE){
    if(length(data) == 1){
      stop("\n More than one species must be supplied \n")
    }
    names_data <- data
  }
  
  utils::data("fishbase")
  list_genus<- fishbase[match(sub("_.*", "", names_data), fishbase$Genus), c("Genus", "Family", "Order")]
  list_local<- data.frame(s= names_data, f= list_genus$Family, o= list_genus$Order)
  not_found<- list_local[which(rowSums(is.na(list_local)) > 0), ]
  
  ####manual insertion##
  print_cat_Family<- function(not_found_fishtree){
    cat("tell the Family of ", not_found_fishtree)
    cat("\n")
  }
  print_cat_Order<- function(not_found_fishtree){
    cat("tell the Order of ", not_found_fishtree)
    cat("\n")
  }
  if(dim(not_found)[1] == 0){
    return(list_local)
  } else{
    for(i in 1:length(not_found$s)){
      spp_family<- readline(prompt = print_cat_Family(not_found_fishtree = not_found$s[i]))
      spp_order<-  readline(prompt = print_cat_Order(not_found_fishtree = not_found$s[i]))
      list_local[which(rowSums(is.na(list_local)) > 0)[i], "f"]<- spp_family
      list_local[which(rowSums(is.na(list_local)) > 0)[i], "o"]<- spp_order
    }
    sub_Perciformes<- which(list_local$o == "Perciformes")
    list_local[sub_Perciformes, "o"]<- "Cichliformes"
    list_local
  }
} 
