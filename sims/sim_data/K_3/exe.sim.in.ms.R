source("../../sim_code/sim.in.ms.R")

K <- 3
if(K == 1){
	ms.command.line <- generate.ms.command.line.values(diploid.population.size = 1e3,
														locus.size = 1e6, 
														per.bp.mu = 1e-8, 
														migration.fraction = 1e-4,
														deep.split = 2e6,
														shallow.split = 1e6)
} else {
	ms.command.line <- generate.ms.command.line.values(diploid.population.size = 1e3,
														locus.size = 1e6, 
														per.bp.mu = 1e-8, 
														migration.fraction = 1e-4,
														deep.split = 2e6,
														shallow.split = 1e6)
}

coords <- as.matrix(expand.grid(1:10,1:10))
N <- nrow(coords)
admix.list <- NULL
if(K > 1){
	w <- sim.admix.props(N,K,coords)
	admix.list <- make.admix.list(N,K,coords,w,0.0001)
}

sim.dataset <- 	generate.conStruct.dataset(n.loci = 1e4,
										K = K,
										coords = coords,
										n.chromo = 20,
										ms.params = ms.command.line,
										admix.list = admix.list,
										subsample = which(coords[,1] > 2 & coords[,1] < 9 &
														   coords[,2] > 2 & coords[,2] < 9))

save(sim.dataset,file="sim.dataset.Robj")


freq.data <- conStruct::process.freq.data(sim.dataset$freq.data$freqs)
D <- fields::rdist(sim.dataset$coords)
pdf(file="ibd.pdf")
	plot(D,freq.data$obsCov)
dev.off()