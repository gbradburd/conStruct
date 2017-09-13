make.ind.genos <- function(x,N){
		geno <- numeric(N)
	if(is.na(x)){
		geno <- rep(-9,N)
	} else {
		geno[sample(1:N,x)] <- 1	
	}
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

load("~/Dropbox/conStruct/data/bears/bear.dataset.redux.Robj")
freqs <- bear.dataset$sample.freqs
sample.sizes <- matrix(2,nrow(freqs),ncol(freqs))
sample.sizes[which(is.na(freqs),arr.ind=TRUE)] <- 0
counts <- freqs*sample.sizes
str.data <- make.structure.data.file(counts,rep(2,nrow(freqs)))
write.table(str.data,file="~/Dropbox/conStruct/data/bears/fastStructure/bears.str",quote=FALSE,row.names=FALSE,col.names=FALSE)
