library(conStruct)
library(doParallel)
library(foreach)
load("sim.dataset.Robj")

options(error=recover)
args <- list("run1" = list("spatial" = FALSE,
						   "geoDist" = fields::rdist(sim.dataset$coords),
						   "K" = 1,
						   "prefix" = "nsp1"),
			 "run2" = list("spatial" = FALSE,
						   "geoDist" = NULL,
						   "K" = 1,
						   "prefix" = "nsp1a"),
			 "run3" = list("spatial" = TRUE,
						   "geoDist" = fields::rdist(sim.dataset$coords),
						   "K" = 1,
						   "prefix" = "sp1"),
			 "run4" = list("spatial" = TRUE,
						   "geoDist" = fields::rdist(sim.dataset$coords),
						   "K" = 3,
						   "prefix" = "sp3"),
			 "run5" = list("spatial" = FALSE,
						   "geoDist" = fields::rdist(sim.dataset$coords),
						   "K" = 3,
						   "prefix" = "nsp3"),
			 "run6" = list("spatial" = FALSE,
						   "geoDist" = NULL,
						   "K" = 3,
						   "prefix" = "nsp3b")
		)

cl <- parallel::makeCluster(4,type="FORK")
doParallel::registerDoParallel(cl)

tmp <- foreach::foreach(i=1:length(args)) %dopar% {
						x <- args[[i]] ; 
        					conStruct::conStruct(spatial = x[["spatial"]],
						  		  			 K = x[["K"]],
						  		  			 freqs = sim.dataset$freq.data$freqs,
						  		  			 geoDist = x[["geoDist"]],
						  		  			 coords = sim.dataset$coords,
						  		  			 prefix = x[["prefix"]])
    			 }

parallel::stopCluster(cl)
