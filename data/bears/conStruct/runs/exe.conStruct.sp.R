library(conStruct)
library(doParallel)
load("../../data/bear.dataset.Robj")

mclapply(1:7,function(k){
	conStruct(spatial = TRUE,
			  K = k,
			  freqs = bear.dataset$sample.freqs,
			  geoDist = fields::rdist.earth(bear.dataset$sample.coords),
			  coords = bear.dataset$sample.coords,
			  prefix = sprintf("bearsK%s_sp",k),
			  n.iter = 5e3)
},mc.cores=7,mc.allow.recursive=FALSE,mc.preschedule=FALSE)