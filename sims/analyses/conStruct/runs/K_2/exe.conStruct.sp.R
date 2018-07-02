library(conStruct)
library(doParallel)
load("../../../../sim_data/K_2/sim.dataset.Robj")

mclapply(1:7,function(k){
	conStruct(spatial = TRUE,
			  K = k,
			  freqs = sim.dataset$freq.data$freqs,
			  geoDist = fields::rdist(sim.dataset$coords),
			  coords = sim.dataset$coords,
			  prefix = sprintf("simK2_K%s_sp",k),
			  n.iter = 5e3)
},mc.cores=6,mc.allow.recursive=FALSE,mc.preschedule=FALSE)