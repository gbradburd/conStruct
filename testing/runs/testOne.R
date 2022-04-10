library(conStruct)

load("sim.dataset.Robj")

options(error=recover)
test <- conStruct::conStruct(spatial = TRUE,
		  		  			 K = 2,
		  		  			 freqs = sim.dataset$freq.data$freqs,
		  		  			 geoDist = fields::rdist(sim.dataset$coords),
		  		  			 coords = sim.dataset$coords,
		  		  			 prefix = "test1")

test <- conStruct::conStruct(spatial = FALSE,
		  		  			 K = 2,
		  		  			 freqs = sim.dataset$freq.data$freqs,
		  		  			 geoDist = fields::rdist(sim.dataset$coords),
		  		  			 coords = sim.dataset$coords,
		  		  			 prefix = "test2",
		  		  			 n.iter=400)

test <- conStruct::conStruct(spatial = FALSE,
		  		  			 K = 2,
		  		  			 freqs = sim.dataset$freq.data$freqs,
		  		  			 geoDist = fields::rdist(sim.dataset$coords),
		  		  			 coords = sim.dataset$coords,
		  		  			 prefix = "test3",
		  		  			 n.iter=500)

test <- conStruct::conStruct(spatial = FALSE,
		  		  			 K = 2,
		  		  			 freqs = sim.dataset$freq.data$freqs,
		  		  			 geoDist = fields::rdist(sim.dataset$coords),
		  		  			 coords = sim.dataset$coords,
		  		  			 prefix = "test4",
		  		  			 n.iter=510)

test <- conStruct::conStruct(spatial = FALSE,
		  		  			 K = 2,
		  		  			 freqs = sim.dataset$freq.data$freqs,
		  		  			 geoDist = fields::rdist(sim.dataset$coords),
		  		  			 coords = sim.dataset$coords,
		  		  			 prefix = "test5",
		  		  			 n.iter=2e3)

test <- conStruct::conStruct(spatial = FALSE,
		  		  			 K = 2,
		  		  			 freqs = sim.dataset$freq.data$freqs,
		  		  			 geoDist = fields::rdist(sim.dataset$coords),
		  		  			 coords = sim.dataset$coords,
		  		  			 prefix = "test5",
		  		  			 n.iter=2e3,
		  		  			 control = setNames(list(0.9),"adapt_delta"))