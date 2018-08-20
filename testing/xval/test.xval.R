library(conStruct)
load("sim.dataset.Robj")

xvals <- x.validation(train.prop = 0.9,
					  n.reps = 4,
					  K = 1:3,
					  freqs = sim.dataset$freq.data$freqs,
					  geoDist = fields::rdist(sim.dataset$coords),
					  coords = sim.dataset$coords,
					  prefix = "xval_test1",
					  n.iter = 1e3,
					  make.figs = FALSE,
					  save.files = FALSE,
					  parallel = TRUE,
					  n.nodes = 4)

save(xvals,file="test1.xvals.Robj")