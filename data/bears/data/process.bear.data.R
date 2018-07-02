################################################################
################################################################
#	script for processing black bear dataset for conStruct analysis
################################################################
################################################################

# .bed, .bim, and .fam file downloaded from dryad repo:
#	http://datadryad.org/resource/doi:10.5061/dryad.dc02b
#
# binary file converted to genotyped counts using:
# 	plink --bfile Puckett_etal_Uamer21k --recodeA --out bear.data --noweb
#
# sampling coordinates are in "Puckett_etal-MaxEntInput" in dryad repo, 
#	but are anonymized w/r/t sample ID, so E. Puckett provided .csv 
#	file associating sample IDs with sampling locations



# format genotype data as allele frequency data matrix
bear.geno <- read.table("bear_data.raw",header=TRUE,stringsAsFactors=FALSE)
sample.names <- bear.geno[,2]
bear.geno <- as.matrix(bear.geno[,-c(1:6)])
bear.freqs <- bear.geno/2

# format metadata
metadata.raw <- read.csv("Puckett-UamerSampleCoordinates.csv",stringsAsFactors=FALSE)
metadata <- as.matrix(metadata.raw)[match(sample.names, metadata.raw$Sample_ID),]
coords <- cbind(as.numeric(metadata[,4]),as.numeric(metadata[,3]))
row.names(coords) <- sample.names

bear.dataset <- list("sample.coords" = coords,
					 "sample.freqs" = bear.freqs,
					 "sample.sizes" = rep(2,nrow(coords)))


# pare down missing data
#	by dropping samples with low genotyping success rates
missing.data <- apply(bear.dataset$sample.freqs,1,function(x){length(which(is.na(x)))})
too.much.missing.data <- which(missing.data/ncol(bear.dataset$sample.freqs) > 0.04)

bear.dataset <- list("sample.freqs" = bear.dataset$sample.freqs[-too.much.missing.data,],
					 "sample.sizes" = bear.dataset$sample.sizes[-too.much.missing.data],
					 "sample.coords" = bear.dataset$sample.coords[-too.much.missing.data,])

save(bear.dataset,file="bear.dataset.Robj")

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

freqs <- bear.dataset$sample.freqs
sample.sizes <- matrix(2,nrow(freqs),ncol(freqs))
sample.sizes[which(is.na(freqs),arr.ind=TRUE)] <- 0
counts <- freqs*sample.sizes
str.data <- make.structure.data.file(counts,rep(2,nrow(freqs)))
write.table(str.data,file="bears.str",quote=FALSE,row.names=FALSE,col.names=FALSE)

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
	genos[which(genos < 0)] <- 0
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

ped_bears <- convert.str.to.ped(str.file="bears.str")
map_bears <- make.map.file(ped_bears)
	write.table(ped_bears,file="../admixture/bears.ped",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(map_bears,file="../admixture/bears.map",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
call <- "plink --noweb --file ../admixture/bears --make-bed --out ../admixture/bears --map3 --allow-no-sex"
system(call)
file.remove("../admixture/bears.ped")
file.remove("../admixture/bears.map")
file.remove("../admixture/bears.log")
file.remove("../admixture/bears.nosex")