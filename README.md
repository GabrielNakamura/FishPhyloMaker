# FishPhyloMaker

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/FishPhyloMaker)](https://cran.r-project.org/package=FishPhyloMaker)

[![codecov](https://codecov.io/gh/GabrielNakamura/FishPhyloMaker/branch/master/graph/badge.svg)](https://codecov.io/gh/GabrielNakamura/FishPhyloMaker)

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)

**Making phylogenies for a pool of local fish species**

The FishPhyloMaker has as the core function `fishPhyloMaker`, that works in a sequential way by download the [most inclusive phylogeny of bony fishes](https://fishtreeoflife.org/) provided by [Rabosky (2018)](https://onlinelibrary.wiley.com/doi/10.1111/jbi.13839) and edit it by replacing and dropping the fish species of the original phylogeny to obtain a phylogeny containing only the desired species specified in a data frame.\
The procedure starts by adding all species of fish that already present any genre in the tree. For those species that do not present any genre, the function finds out for all species of the same family of the species that must be inserted. By an interactively procedure, the user must specify which species from that family in the tree the species to be inserted is most related. Finally, for those species that do not present any family representatives in the tree their inclusion is made by searching for the species of the same order presented in the original tree. If any, the functions asks the user to provide a newick file to be binded to the tree or to insert the species as politomies in the order of species. The function returns a newick file containing the phylogeny of species provide in data argument.

The user must provide to `fishPhyloMaker` function a data frame that present the following format:

| genre  | species |  family |  order |
|--------|:-------:|:-------:|:------:|
| Genre1 |   sp1   | Family1 | Order1 |
| Genre2 |   sp2   | Family2 | Order2 |
| Genre3 |   sp3   | Family3 | Order3 |


This table can be done mannually or by passing to `tab_function` a list of species or a community data matrix with species names in columns.


```
