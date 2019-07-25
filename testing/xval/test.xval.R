library(conStruct)
load("sim.dataset.Robj")

options(error=recover)

xvals <- x.validation(train.prop = 0.9,
					  n.reps = 2,
					  K = 1:2,
					  freqs = sim.dataset$freqs,
					  geoDist = fields::rdist(sim.dataset$coords),
					  coords = sim.dataset$coords,
					  prefix = "xval_test1",
					  n.iter = 1e3,
					  make.figs = FALSE,
					  save.files = FALSE,
					  parallel = FALSE,
					  n.nodes = NULL)

save(xvals,file="test1.xvals.Robj")