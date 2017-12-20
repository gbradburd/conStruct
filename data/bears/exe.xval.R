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
					},mc.cores=6,mc.allow.recursive=FALSE,mc.preschedule=FALSE)
	names(x.val) <- paste0("rep_",1:n.reps)
	return(x.val)
}


load("bear.dataset.Robj")

missing.data <- apply(bear.dataset$sample.freqs,1,function(x){length(which(is.na(x)))})
too.much.missing.data <- which(missing.data/ncol(bear.dataset$sample.freqs) > 0.04)

bear.dataset <- list("sample.freqs" = bear.dataset$sample.freqs[-too.much.missing.data,],
					 "sample.sizes" = bear.dataset$sample.sizes[-too.much.missing.data],
					 "sample.coords" = bear.dataset$sample.coords[-too.much.missing.data,])

save(bear.dataset,file="bear.dataset.Robj")
x.validation(test.pct = 0.10,
			 n.reps = 10,
			 K = 1:7,
			 freqs = bear.dataset$sample.freqs,
			 geoDist = fields::rdist.earth(bear.dataset$sample.coords),
			 coords = bear.dataset$sample.coords,
			 prefix = "bear_",
			 n.iter = 5e3)
