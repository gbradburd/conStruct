library(conStruct)
load("sim.dataset.Robj")

library(foreach)
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

xvals <- x.validation(train.prop = 0.9,
					  n.reps = 4,
					  K = 1:4,
					  freqs = sim.dataset$freq.data$freqs,
					  data.partitions = NULL,
					  geoDist = fields::rdist(sim.dataset$coords),
					  coords = sim.dataset$coords,
					  prefix = "xval_test2",
					  n.iter = 1e3,
					  make.figs = FALSE,
					  save.files = FALSE,
					  parallel = TRUE,
					  n.nodes = 4)

save(xvals,file="xvals2.Robj")

stopCluster(cl)