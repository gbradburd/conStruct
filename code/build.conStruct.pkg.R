setwd("~/Dropbox/conStruct/code/conStruct")
require(devtools)
require(rmarkdown)
file.remove(list.files(pattern=".tar.gz"))

#setup() #ONLY RUN setup() once!
# edit description file
# put R code in R file
# add dependencies
oneK.stan.block <- paste(readLines("~/Dropbox/conStruct/code/stan_model_blocks/basic.txt",warn=FALSE),collapse="\n")
space.oneK.stan.block <- paste(readLines("~/Dropbox/conStruct/code/stan_model_blocks/space.txt",warn=FALSE),collapse="\n")
multiK.stan.block <- paste(readLines("~/Dropbox/conStruct/code/stan_model_blocks/multiK.txt",warn=FALSE),collapse="\n")
space.multiK.stan.block <- paste(readLines("~/Dropbox/conStruct/code/stan_model_blocks/space_multiK.txt",warn=FALSE),collapse="\n")
use_data(oneK.stan.block,
		 space.oneK.stan.block,
		 multiK.stan.block,
		 space.multiK.stan.block,
		 internal=TRUE,overwrite=TRUE)


# edit NAMESPACE
load_all()
# write documentation
document()
# after updating documentation, run document() again, then run
# R CMD Rdconv -t html manfile.Rd -o out.html
#have to add texbin to R's path so that it can find pdflatex
#	to build the vignette
Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/usr/texbin",sep=":"))
check(build_args="--resave-data")
build(path="~/Dropbox/conStruct/code/conStruct",args="--resave-data")
install()