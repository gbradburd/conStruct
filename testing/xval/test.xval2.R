library(conStruct)
load("sim.dataset.Robj")

library(foreach)
library(doParallel)
cl <- makeCluster(2,type="FORK")
registerDoParallel(cl)

xvals <- x.validation(train.prop = 0.9,
					  n.reps = 2,
					  K = 1:2,
					  freqs = sim.dataset$freqs,
					  data.partitions = NULL,
					  geoDist = fields::rdist(sim.dataset$coords),
					  coords = sim.dataset$coords,
					  prefix = "xval_test2",
					  n.iter = 1e3,
					  make.figs = FALSE,
					  save.files = FALSE,
					  parallel = TRUE,
					  n.nodes = 2)

save(xvals,file="xvals2.Robj")

stopCluster(cl)