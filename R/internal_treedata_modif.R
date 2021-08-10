#' Internal function to show the number of species from species pool that lacks in phylogeny
#'
#' @param phy Phylogenetic hypothesis in newick format
#' @param data Data frame with species to be added in tree
#' @param sort Sorting species in alphabetic order
#' @param warnings Logical
#'
#' @return character vector
#' @keywords internal
#' 
treedata_modif<- function (phy, data, sort = FALSE, warnings = TRUE) 
{
  dm = length(dim(data))
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  if (is.factor(data)) {
    data <- as.matrix(data)
  }
  if (is.array(data) & length(dim(data)) == 1) {
    data <- as.matrix(data)
  }
  if (is.null(rownames(data))) {
    stop("names for 'data' must be supplied")
  }
  else {
    data.names <- rownames(data)
  }
  nc <- geiger::name.check(phy, data)
  if (is.na(nc[[1]][1]) | nc[[1]][1] != "OK") {
    if (length(nc[[1]] != 0)) {
      phy = ape::drop.tip(phy, as.character(nc[[1]]))
      if (warnings) {
        warning(paste("The following tips were not found in 'data' and were dropped from 'phy':\n\t", 
                      paste(nc[[1]], collapse = "\n\t"), sep = ""))
      }
    }
    if (length(nc[[2]] != 0)) {
      m <- match(data.names, nc[[2]])
      data = as.matrix(data[is.na(m), ])
      data.names <- data.names[is.na(m)]
      if (warnings) {
        warning(paste("The following tips were not found in 'phy' and were dropped from 'data':\n\t", 
                      paste(nc[[2]], collapse = "\n\t"), sep = ""))
      }
    }
  }
  order <- match(data.names, phy$tip.label)
  rownames(data) <- phy$tip.label[order]
  if (sort) {
    index <- match(phy$tip.label, rownames(data))
    data <- as.matrix(data[index, ])
  }
  if (dm == 2) {
    data <- as.matrix(data)
  }
  phy$node.label = NULL
  return(list(phy = phy, data = data, nc= nc))
}