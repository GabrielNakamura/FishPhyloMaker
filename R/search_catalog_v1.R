#' @title Searching `Eschmeyer's Catalog of Fish`
#' @description This function get fields for fetching valid species information from the `Eschmeyer's Catalog of Fishes`.
#' @rdname search_cas
#' @name search_cas
#' @param query \string One and more queries, support family, genus, and species.
#' @param type \string There is one required parameter, type = c("genus_family","species_family","species_genus","species").For example, `genus_family` for genera in a family,`species_family` for species in a family, `species_genus` for species in a genus, and `species` for species related to a species.
#' @format The following data frame of the fields to be returned are recognized:
#' \itemize{
#' \item query : family, genus, and species Matched.
#' \item species(genus)_author : species(genus) and author(s) returned.
#' \item species(genus) : Species(genus) returned.
#' \item author : Author(s) returned.
#' \item family : Family(subfamily) of species.
#' \item status : Current status including Validation, synonym, and Uncertainty.
#' }
#' @return Data frame.
#' @author Liuyong Ding \email{ly_ding@126.com}
#' @importFrom httr GET
#' @importFrom httr status_code
#' @importFrom magrittr `%>%`
#' @importFrom xml2 read_html()
#' @importFrom xml2 xml_find_all
#' @importFrom xml2 xml_text
#' @importFrom tibble enframe()
#' @importFrom dplyr mutate
#' @importFrom stringr str_squish
#' @importFrom dplyr filter
#' @importFrom stringr str_extract
#' @importFrom stringr str_split
#' @importFrom stringr str_trim
#' @importFrom tibble as_tibble
#' @details See for website \url{http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp} details.
#' @examples
#' \dontrun{
#' # For genera in a family:
#' r1 <- search_cas(query = c("cyprinidae","balitoridae"),type = "genus_family")
#'
#' # For species in a family:
#' r2 <- search_cas(query = "balitoridae",type = "species_family")
#'
#' # For species in a genus:
#' r3 <- search_cas(query = "Brachyplatystoma",type = "species_genus")
#' 
#' # For species related to a species or subspecies:
#' r4 <- search_cas(query = c("Hoplias malabaricus", "Brachyplatystoma filamentosum"),type = "species")
#' 
#' }
#' @export
#' 


search_cas <- function(query,type){
  # multiple species query 
  if(length(query) > 1){
    res <- lapply(query, function(x){
      get_cas(x,type = type)
    })
    res <- do.call("rbind",res)
    
  }else{
    res <- get_cas(query,type = type)
  }
  return(res)
}

get_cas <- function(query,type){
  baseurl = "https://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatget.asp?"
  stopifnot(is.character(query))
  #type = c("genus_family","species_family","species_genus","species")
  switch(type, `genus_family` = {
    url = paste0(baseurl,"tbl=genus&family=",query)
  }, `species_family` = {
    url = paste0(baseurl,"tbl=species&family=",query)
  }, `species_genus` = {
    url = paste0(baseurl,"tbl=species&genus=",query)
  },`species` = {
    query = gsub(" ","_",query)
    url = paste0(baseurl,"tbl=species&genus=",str_split(query,"_",simplify = T)[1],"&species=",str_split(query,"_",simplify = T)[2])
  }
  )
  res = httr::GET(url)
  if (httr::status_code(res) == 200){
    result = res %>%
      xml2::read_html() %>%
      xml2::xml_find_all("//p[@class='result']") %>%
      xml2::xml_text() %>%
      tibble::enframe() %>%
      dplyr::mutate(value = stringr::str_squish(value)) %>%
      dplyr::filter(value != "") %>%
      dplyr::mutate(
        query = stringr::str_extract(value, "[a-z]+, [A-Za-z]+"),
        content = stringr::str_extract(value, "(?<=Current status: ).*"),
        species = stringr::str_split(content,"\\. ",simplify = T)[,1],
        family = stringr::str_split(content,"\\. ",simplify = T)[,2],
        species_author = stringr::str_trim(ifelse(grepl("Valid as",species),gsub("Valid as ","",species),
                                                  ifelse(grepl("Synonym of",species),gsub("Synonym of ","",species),
                                                         ifelse(grepl("Uncertain as",species),gsub("Uncertain as ","",species),species)))),
        status = ifelse(grepl("Valid as",species),"Validation",
                        ifelse(grepl("Synonym of",species),"Synonym",
                               ifelse(grepl("Uncertain as",species),"Uncertainty",NA)))
      ) %>%
      na.omit() %>%
      dplyr::select(query,species_author,family,status)
    
    dd = str_locate_all(result$species_author, " ")
    end = c()
    switch(type, `genus_family` = {
      for (i in 1:length(dd)) {
        end[i] = dd[[i]][1]
      }
    }, `species_family` = {
      for (i in 1:length(dd)) {
        end[i] = dd[[i]][2]
      }
    }, `species_genus` = {
      for (i in 1:length(dd)) {
        end[i] = dd[[i]][2]
      }
    },`species` = {
      for (i in 1:length(dd)) {
        end[i] = dd[[i]][2]
      }
    }
    )
    result$species = str_sub(result$species_author, 1, end - 1)
    result$author = str_sub(result$species_author, end + 1)
    if(type == "genus_family"){
      names(result)[2] = "genus_author"
      names(result)[5] = "genus"
      result$family = gsub("\\.","",result$family)
    }
    return(tibble::as_tibble(result[,c(1,2,5,6,3,4)]))
  }else{
    cat("Error request - the parameter query is not valid")
    browseURL("http://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp")
  }
  
}
