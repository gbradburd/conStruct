library(conStruct)
library(doParallel)
library(foreach)
load("sim.dataset.Robj")

args <- list("run1" = list("spatial" = FALSE,
						   "geoDist" = fields::rdist(sim.dataset$coords),
						   "K" = 1,
						   "prefix" = "nsp1"),
			 "run2" = list("spatial" = TRUE,
						   "geoDist" = fields::rdist(sim.dataset$coords),
						   "K" = 1,
						   "prefix" = "sp1"),
			 "run3" = list("spatial" = FALSE,
						   "geoDist" = fields::rdist(sim.dataset$coords),
						   "K" = 3,
						   "prefix" = "nsp3"),
			 "run4" = list("spatial" = TRUE,
						   "geoDist" = fields::rdist(sim.dataset$coords),
						   "K" = 3,
						   "prefix" = "sp3")
		)

cl <- parallel::makeCluster(4)
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
