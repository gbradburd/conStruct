make.ind.genos <- function(ones,N){
	geno <- numeric(N)
	geno[sample(1:N,ones)] <- 1
	return(geno)
}

split.pop.into.inds <- function(pop.counts,pop.nChromo){
	ind.geno <- apply(pop.counts,2,function(x){make.ind.genos(x,pop.nChromo)})
	return(ind.geno)
}

make.structure.data.file <- function(counts,sample.sizes){
	geno.matrix <- matrix(NA,sum(sample.sizes),ncol(freqs))
	index <- 0
	for(n in 1:nrow(freqs)){
		geno.matrix[index+c(1:sample.sizes[n]),] <- split.pop.into.inds(counts[n,,drop=FALSE],sample.sizes[n])
		index <- index + sample.sizes[n]
	}
	PopData <- unlist(lapply(1:length(sample.sizes),
					function(n){
						lapply(1:sample.sizes[n],
							function(i){
								sprintf("ind_%s.%s",n,i)
							})
					}))
	LocData <- unlist(lapply(1:length(sample.sizes),function(n){rep(n,sample.sizes[n])}))
	return(cbind(PopData,LocData,geno.matrix))
}

load("~/Dropbox/conStruct/sims/cross_validation/K_1/sim.dataset.Robj")
freqs <- sim.dataset$freq.data$freqs
sample.sizes <- sim.dataset$freq.data$sample.sizes
counts <- freqs*sample.sizes
str.data <- make.structure.data.file(counts,sample.sizes[,1])
write.table(str.data,file="datasets/simK1/simK1.str",quote=FALSE,row.names=FALSE,col.names=FALSE)

load("~/Dropbox/conStruct/sims/cross_validation/K_2/sim.dataset.Robj")
freqs <- sim.dataset$freq.data$freqs
sample.sizes <- sim.dataset$freq.data$sample.sizes
counts <- freqs*sample.sizes
str.data <- make.structure.data.file(counts,sample.sizes[,1])
write.table(str.data,file="datasets/simK2/simK2.str",quote=FALSE,row.names=FALSE,col.names=FALSE)

load("~/Dropbox/conStruct/sims/cross_validation/K_3/sim.dataset.Robj")
freqs <- sim.dataset$freq.data$freqs
sample.sizes <- sim.dataset$freq.data$sample.sizes
counts <- freqs*sample.sizes
str.data <- make.structure.data.file(counts,sample.sizes[,1])
write.table(str.data,file="datasets/simK3/simK3.str",quote=FALSE,row.names=FALSE,col.names=FALSE)