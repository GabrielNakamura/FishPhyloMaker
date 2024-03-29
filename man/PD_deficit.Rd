% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PD_deficit.R
\name{PD_deficit}
\alias{PD_deficit}
\title{Title Calculate the amount of phylogenetic deficit in assemblages}
\usage{
PD_deficit(phylo, data, level = "Congeneric_insertion")
}
\arguments{
\item{phylo}{Phylogenetic tree in newick format, can be an object from \code{\link{FishPhyloMaker}} function}

\item{data}{A data frame containing the classification  informing the level of insertions. This can be obtained
from \code{\link{FishPhyloMaker}} function}

\item{level}{Character indicating which level must be considered in the calculation of PD deficit.
Can be a vector with the levels ("Congeneric_insertion", "Congeneric_Family_level", "Family_insertion", "Order_insertion")
which will be considered in the calculation of phylogenetic deficit.
default is "Congeneric_insertion".}
}
\value{
A vector containing four values:

\if{html}{\out{<div class="sourceCode">}}\preformatted{- Amount phylogenetic information present in the tree before insertions (PDintree)

- Amount of phylogenetic information inserted in the tree (PDdeficit)

- Total Phylogenetic information of the tree (PDtotal)

- A ratio calculated as PDdeficit/PDtotal (Darwinian_deficit)
}\if{html}{\out{</div>}}
}
\description{
Title Calculate the amount of phylogenetic deficit in assemblages
}
\seealso{
\code{\link{FishPhyloMaker}} for phylogeny and data frame containing the classification of insertions
}
