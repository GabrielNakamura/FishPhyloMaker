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

MANUAL CHECKING

* Please omit the redundant "An R Package to" from your title.
  
  I removed the redundant part of the title
  
* Is there a doi available to the reference in your description field you could add in the form <doi:10.prefix/suffix>?
  
  I added a doi to the reference
  
* Please unwrap the examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing.
(You could also replace \dontrun{} with \donttest, if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions. Otherwise, you can also write some tests.)
  
  I remove  \dontrun{} and change it to \donttest{} in fishPhyloMaker.R and whichFishAdd.R functions, since the example are executed in > 5 sec. I can not remove \dontrun{} from tab_function.R since it depends on a interactive procedure with the user.  


* Possibly mis-spelled words in DESCRIPTION:
  al (19:91)
  et (19:88)
  Phylogenies (3:8)
  Rabosky (19:80)
  
The spelling is correct for all words


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