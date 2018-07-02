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
poplar.sample.sizes <- numeric(length(drainages))

# function for getting mean number of samples genotyped 
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
	poplar.sample.sizes[i] <- length(in.pop)
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
save(poplar.sample.sizes,file="poplar.sample.sizes.Robj")

################################
# make STRUCTURE and ADMIXTURE input data files
################################

################################
#	structure
################################
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

freqs <- poplar.data$freqs
sample.sizes <- poplar.sample.sizes*2
counts <- freqs*matrix(sample.sizes,nrow(freqs),ncol(freqs))
str.data <- make.structure.data.file(counts,sample.sizes)
write.table(str.data,file="poplars.str",quote=FALSE,row.names=FALSE,col.names=FALSE)


################################
#	admixture
################################
convert.str.to.ped <- function(str.file){
	str <- read.table(str.file,header=FALSE,stringsAsFactors=FALSE,row.names=1)
	n.loci <- ncol(str)-1
	n.haps <- nrow(str)
	hap.IDs <- row.names(str)
	pop.IDs <- as.numeric(unlist(lapply(lapply(strsplit(hap.IDs,"_"),"[[",2),function(x){strsplit(x,"\\.")[[1]][[1]]}))[seq(1,n.haps,by=2)])
	mand.cols <- matrix(0,nrow=n.haps/2,ncol=6)
	#family
	mand.cols[,1] <- pop.IDs
	#ind ID
	mand.cols[,2] <- unlist(lapply(unique(pop.IDs),function(x){1:length(which(pop.IDs==x))}))
	genos <- matrix(NA,nrow=n.haps/2,ncol=n.loci*2)
	genos[,seq(1,n.loci*2,by=2)] <- as.matrix(str[seq(1,n.haps,by=2),2:ncol(str)])
	genos[,seq(2,n.loci*2,by=2)] <- as.matrix(str[seq(2,n.haps,by=2),2:ncol(str)])
	genos <- genos + 1
	ped <- cbind(mand.cols, genos)
	return(ped)
}

make.map.file <- function(ped){
	#recover()
	n.snps <- (ncol(ped)-6)/2
	map <- matrix(NA,nrow=n.snps,ncol=3)
	#CHR
	map[,1] <- rep(1,n.snps)
	#RS
	map[,2] <- paste0("rs",1:n.snps)
	#BP
	map[,3] <- seq(1,n.snps*1e3,length.out=n.snps)
	return(map)
}

ped_poplar <- convert.str.to.ped(str.file="poplars.str")
map_poplar <- make.map.file(ped_poplar)
	write.table(ped_poplar,file="../admixture/poplar.ped",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(map_poplar,file="../admixture/poplar.map",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
call <- "plink --noweb --file ../admixture/poplar --make-bed --out ../admixture/poplar --map3 --allow-no-sex"
system(call)
file.remove("../admixture/poplar.ped")
file.remove("../admixture/poplar.map")
file.remove("../admixture/poplar.log")
file.remove("../admixture/poplar.nosex")