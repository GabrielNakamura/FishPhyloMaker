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
   URL: %22http://www.fishbase.org%22
     From: inst/doc/FishPhyloMaker_vignette.html
     Message: Invalid URI scheme
     
I changed the URL to [Fishbase database](https://www.fishbase.se) and now the link is correct with the proper https address

* Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1590/1982-0224-2020-0126
    From: man/neotropical_comm.Rd
    Status: Error
    Message: libcurl error code 60:
      	SSL certificate problem: unable to get local issuer certificate
      	(Status without verification: OK)

The URL is valid 