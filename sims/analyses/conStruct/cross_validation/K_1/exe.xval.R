library(conStruct)
library(doParallel)

x.validation <- function(test.pct,n.reps,K,freqs,geoDist,coords,prefix,n.iter){
	x.val <- mclapply(1:n.reps,
					function(i){
						x.validation.rep(rep.no = i,
										 test.pct,
										 K,
										 freqs,
										 geoDist,
										 coords,
										 prefix,
										 n.iter)
					},mc.cores=5)
	names(x.val) <- paste0("rep_",1:n.reps)
	return(x.val)
}


load("../../../../sim_data/K_1/sim.dataset.Robj")

x.validation(test.pct = 0.10,
			 n.reps = 10,
			 K = 1:7,
			 freqs = sim.dataset$freq.data$freqs,
			 geoDist = fields::rdist.earth(sim.dataset$coords),
			 coords = sim.dataset$coords,
			 prefix = "simK1_",
			 n.iter = 1e3)