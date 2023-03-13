
## conStruct (continuous Structure) ReadMe

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/gbradburd/conStruct/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gbradburd/conStruct/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->


This repo contains the code for the method **conStruct** - a statistical tool
for modeling continuous and discrete population genetic structure.

The manuscript, data files, and analysis scripts associated with the publication
"Inferring Continuous and Discrete Population Genetic Structure Across Space,"
have been moved, and can be accessed at the links below:

 * [paper](https://doi.org/10.1534/genetics.118.301333)
 * [manuscript repo](https://github.com/gbradburd/conStruct-paper)
 * [data dryad](https://doi.org/10.5061/dryad.5qj7h09)

## Installation

### Latest release

To install the most recent release of the **conStruct** R package:

```r
install.packages("conStruct")
```

Upon installation, the **conStruct** models will be compiled, which may
spit lots of text, and possibly some warnings, to your screen. This is
totally normal, and you should only be concerned if you get errors
and the installation fails.


### Development version

To install the development version from github:

```r
	library(devtools)
	install_github("gbradburd/conStruct",build_vignettes=TRUE)
```

Note that Windows users may have to download Rtools as a
standalone executable before trying to install the **conStruct** R package.


## Getting Started

A complete manual for all documented functions is available [here](https://github.com/gbradburd/conStruct/blob/master/man/conStruct-manual.pdf).

In addition, there are four vignettes included in the package that walk through
various steps in the analysis pipeline in detail. You can find them using:

```r
# formatting data
	vignette(topic="format-data",package="conStruct")

# how to run a conStruct analysis
	vignette(topic="run-conStruct",package="conStruct")

# how to visualize the output of a conStruct model
	vignette(topic="visualize-results",package="conStruct")

# how to compare and select between different conStruct models
	vignette(topic="model-comparison",package="conStruct")
```

There is also an example data file included in the package, which you can
load using the command:

```r
	data(conStruct.data)
```

## Contact

After referring to the manual and vignettes,
please direct all queries to bradburd (at) umich.edu,
or post as issues on the git repo.