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