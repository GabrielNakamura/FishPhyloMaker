## Test environments

* local OS X install, R 4.1.0
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Fedora Linux, R-devel, clang, gfortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC

## R CMD check results

There were no ERRORs or WARNINGs.

NOTES.

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Gabriel Nakamura <gabriel.nakamura.souza@gmail.com>’

New submission

* Possibly misspelled words in DESCRIPTION:
  Rabosky (20:80)
  
The spelling is correct

* Found the following (possibly) invalid URLs:
  URL: https://www.scielo.br/scielo.php?script=sci_arttext&pid=S1679-62252021000100211 (moved to https://www.scielo.br/j/ni/a/jWsfj3ZDr6qBbMyg5kfczrB/?lang=en)
From: man/neotropical_comm.Rd

The URL is valid
    
* Found the following (possibly) invalid URLs:
   URL: %22https://www.fishbase.se/search.php%22
     From: inst/doc/FishPhyloMaker_vignette.html
     Message: Invalid URI scheme
     
I changed the URL to https://www.fishbase.org

* URL: https://cran.r-project.org/web/packages/rfishbase/rfishbase.pdf
     From: README.md
     Status: 200
     Message: OK
     CRAN URL not in canonical form
   The canonical URL of the CRAN page for a package is
     https://CRAN.R-project.org/package=pkgname
     
I changed the URL to the canonical URL of the CRAN page https://CRAN.R-project.org/package=rfishbase

*  Found the following URLs which should use \doi (with the DOI name only):
   File 'neotropical_comm.Rd':
     https://doi.org/10.1590/1982-0224-2020-0126
    
I changed the URL to \doi{10.1590/1982-0224-2020-0126}

