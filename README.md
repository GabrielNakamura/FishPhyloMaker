# FishPhyloMaker

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/FishPhyloMaker)](https://cran.r-project.org/package=FishPhyloMaker)

[![codecov](https://codecov.io/gh/GabrielNakamura/FishPhyloMaker/branch/master/graph/badge.svg)](https://codecov.io/gh/GabrielNakamura/FishPhyloMaker)

[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://choosealicense.com/licenses/mit/)

**Making phylogenies for a pool of local fish species**

The FishPhyloMaker package has as the core function `FishPhyloMaker`, that works by downloading the [most inclusive phylogeny of bony fishes](https://fishtreeoflife.org/) provided by [Rabosky (2018)](https://onlinelibrary.wiley.com/doi/10.1111/jbi.13839) and edit this phylogeny, in a sequential way by replacing and dropping the fish species of the original phylogeny to obtain a phylogenetic tree containing only the desired species specified in a data frame. The procedure starts by adding all species of fish that already present any Genus in the tree (Congeneric species). For species that do not present any representative of the same Genus, the function finds out for all species of the same family, if not present any species in the tree the search is performed for species of the same Order. By an interactively procedure, the user must specify which Genus from that family (or families in orders) in the tree the species to be inserted is most related. The function returns a newick file containing the phylogeny of species provide in data argument and a data frame (if argument `return.insertions = TRUE`) containing all species with an character indicating at which level the species was added in the tree. There are four possible categories in which species can be tagged regarding their order of insertion:

-   **Present_in_tree** the species was already present in the original tree;
-   **Congeneric_insertion** species inserted as a sister species of the same genus presented in the tree;
-   **Family_insertion** species inserted as being related to a given genus, between two genus or at the node that corresponds to the Family;
-   **Order_insertion** species inserted as being related to a given family, between two families or at the node that corresponds to the Family;

The user must provide to `FishPhyloMaker` function a data frame that present the following format:

|    s  |    f    |    o   |
|:-----:|:-------:|:------:|
| G_sp1 | Family1 | Order1 |
| G_sp2 | Family2 | Order2 |
| G_sp3 | Family3 | Order3 |

This table can be done mannually or by passing to `tab_function` a list of species or a community data matrix with species names in columns.

To install the package the user must type:

```{r downpkg, eval=F, echo = F}

devtools::install_github("GabrielNakamura/FishPhyloMaker", ref = "main")

```

To run an example the user can load a dataset contained in the package:

```{r examp, eval=T, echo = F}

data(neotropical_comm)
data_comm <- neotropical_comm[, -c(1, 2)] # removing latitude and longitude


```

First the user must obtain the data necessary to enter in `FishPhyloMaker` function using `tab_function`

```{r tab_examp, eval=T, echo = T}

taxon_data <- tab_function(data_comm)

```

And finally run `FishPhyloMaker`

```{r maker_examp, eval=T, echo = T}

res_phylo <- FishPhyloMaker(data = taxon_data, return.insertions = TRUE)

```

The output has two objects, a newick that contains the phylogeny for the local pool of species, and can be
   directly plot

```{r plot_examp, eval=T, echo = T}

plot(res_phylo$Phylogeny, cex = 0.7)

```
And a data frame containing in which level the species of local pool was inserted

```{r table_examp, eval=T, echo = T}

res_phylo$Insertions_data

```
