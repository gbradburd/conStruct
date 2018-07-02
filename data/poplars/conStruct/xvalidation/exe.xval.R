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
					},mc.cores=7,mc.allow.recursive=FALSE,mc.preschedule=FALSE)
	names(x.val) <- paste0("rep_",1:n.reps)
	return(x.val)
}


load("../../data/poplar.data.Robj")

x.validation(test.pct = 0.10,
			 n.reps = 10,
			 K = 1:7,
			 freqs = freqs,
			 geoDist = fields::rdist.earth(poplar.data$coords),
			 coords = poplar.data$coords,
			 prefix = "poplar_",
			 n.iter = 2e3)