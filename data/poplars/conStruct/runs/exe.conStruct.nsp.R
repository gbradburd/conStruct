library(conStruct)
library(doParallel)
load("../../data/poplar.data.Robj")

mclapply(1:7,function(k){
	conStruct(spatial = FALSE,
			  K = k,
			  freqs = poplar.data$freqs,
			  geoDist = fields::rdist.earth(poplar.data$coords),
			  coords = poplar.data$coords,
			  prefix = sprintf("poplarsK%s_nsp",k),
			  n.iter = 5e3)
},mc.cores=6,mc.allow.recursive=FALSE,mc.preschedule=FALSE)

