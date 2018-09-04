
## Test environments
* local OS X install R 3.5.0
* win-builder

## R CMD check results
There were no ERRORs or WARNINGs

There were 2 NOTES:

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.

Explanation: GNU make is required for packages that 
depend on rstan and which are developed using rstantools. 
The requirement is noted in the DESCRIPTION file.

* checking installed package size ... NOTE
  installed size is 11.7Mb
  sub-directories of 1Mb or more:
    data   6.6Mb
    libs   4.1Mb

Explanation: This package has a large installed library 
size because it uses the Stan MCMC library as a backend, 
and the C++ Stan models included in the package are 
compiled upon installation. It has a large data directory 
because it includes three example files that are used in 
the vignettes: a file that contains example data inputs 
to run an analysis, and two example output files from an 
analysis, which are used to demonstrate post-analysis 
visualization functions.


## Downstream dependendencies

* There are currently no downstream dependencies for this package.