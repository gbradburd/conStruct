################################################################
################################################################
#	script for processing poplar dataset for conStruct analysis
################################################################
################################################################

##downloaded from http://datadryad.org/resource/doi:10.5061/dryad.7s848
poplar <- read.table("FOR_DRYAD_TandB_fwd_june2013_434.txt",as.is=TRUE,sep="\t")

## from suppmat tables from the paper
pop.metadata <- read.csv(file="Populus_metadata.csv",as.is=TRUE,skip=1)

# get samples for which there's genetic data
# 	only retain metadata for genotyped samples
pop.samples <- poplar[10,]
map <- match(pop.samples,pop.metadata$Accession)
pop.metadata <- pop.metadata[map,]
pop.metadata <- pop.metadata[-1,]

# convert genotype nucleotide data into 
#	allele frequency data matrix for all individuals
poplar.genos <- poplar[11:nrow(poplar),2:ncol(poplar)]
poplar.individual.freqs <- apply(poplar.genos,1,
									function(geno){
										genotype <- sapply(as.character(geno),function(x){strsplit(x,"/|")[[1]][1:2]})
										bases <- unique(c(genotype[1,],genotype[2,]))
										my.bases <- bases[bases %in% c("A","T","C","G")]
										my.bases <- sample(my.bases)
										allele.1 <- (genotype[1,] == my.bases[1]) + (genotype[2,] == my.bases[1]) 
										allele.2 <- (genotype[1,] == my.bases[2]) + (genotype[2,] == my.bases[2]) 
										return(allele.1 / (allele.1+allele.2))
								})

# collapse individual data into drainage-specific populations
drainages <- unique(pop.metadata$Drainage.Location.name)
drainage.freqs <- matrix(NA,nrow=length(drainages),ncol=ncol(poplar.individual.freqs))
drainage.sample.sizes <- matrix(NA,nrow=length(drainages),ncol=ncol(poplar.individual.freqs))
drainage.long <- numeric(length(drainages))
drainage.lat <- numeric(length(drainages))
drainage.sp <- numeric(length(drainages))

# function for getting number of samples genotyped 
#	in each drainage at each locus
get.pop.sample.sizes <- function(ind.freqs,in.pop){
	pop.sample.sizes <- rep(length(in.pop),ncol(ind.freqs)) - 
							apply(ind.freqs[in.pop,,drop=FALSE],2,
									function(x){
										length(which(is.na(x)))
									})
	return(2 * pop.sample.sizes)
}

# loop through drainages and average over individual
#	genotype and sample location
for(i in 1:length(drainages)){
	in.pop <- which(pop.metadata$Drainage.Location.name==drainages[i])
	drainage.freqs[i,] <- colMeans(poplar.individual.freqs[in.pop,,drop=FALSE],na.rm=TRUE)
	drainage.sample.sizes[i,] <- get.pop.sample.sizes(poplar.individual.freqs,in.pop)
	drainage.long[i] <- mean(as.numeric(pop.metadata$Longitude[in.pop]))
	drainage.lat[i] <- mean(as.numeric(pop.metadata$Latitude[in.pop]))
	if(length(unique(pop.metadata$Species[in.pop])) < 2){
		drainage.sp[i] <- unique(pop.metadata$Species[in.pop])
	}	else {
		drainage.sp[i] <- "mixed"
	}
}

drainage.coords <- cbind(drainage.long,drainage.lat)
rownames(drainage.coords) <- unique(pop.metadata$Drainage.Location.name)


# drop all missing data
missing.data <- apply(drainage.freqs,2,function(x){any(is.na(x))})
drainage.freqs <- drainage.freqs[,!missing.data]
drainage.sample.sizes <- drainage.sample.sizes[,!missing.data]


poplar.data <- list("freqs" = drainage.freqs,
					"coords" = drainage.coords,
					"sample.sizes" = rowMeans(drainage.sample.sizes,na.rm=TRUE),
					"sp.ID" = drainage.sp)

save(poplar.data,file="poplar.data.Robj")
