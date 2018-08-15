#' Convert a dataset from STRUCTURE to conStruct format
#'
#' \code{structure2conStruct} converts a STRUCTURE dataset 
#' to conStruct format
#' 
#' This function takes a population genetics dataset in 
#' STRUCTURE format and converts it to conStruct format. 
#' The STRUCTURE file must have one row per individual 
#' and two columns per locus, and can only contain 
#' bi-allelic SNPs.
#' 
#' @param infile The name and path of the file in STRUCTURE format 
#' 			to be converted to \code{conStruct} format. 
#' @param start.loci The index of the first column in the dataset 
#'			that contains genotype data.
#' @param missing.datum The character or value used to denote 
#' 			missing data in the STRUCTURE dataset (often -9).
#' @param outfile The name and path of the file containing the 
#'			\code{conStruct} formatted dataset to be generated 
#' 			by this function.
#'
#' @details This function takes a STRUCTURE format data file and 
#'		converts it to a \code{conStruct} format data file. 
#' 		The STRUCTURE dataset should be in the ONEROWPERIND 
#' 		file format, with one row per individual and two columns 
#' 		per locus (this function therefore can only be applied to 
#' 		diploid organisms). The first column of the STRUCTURE dataset 
#' 		should be individual names. There may be any number of other 
#' 		columns that contain non-genotype information before the first
#'		column that contains genotype data, but there can 
#' 		be no extraneous columns at the end of the dataset, after the 
#' 		genotype data.  The genotype data should be bi-allelic 
#'		single nucleotide polymorphisms (SNPs).
#'		
#'	@return This function returns an allele frequency data matrix 
#'		that can be used as the \code{freqs} argument in a conStruct 
#'		analysis run using \code{\link{conStruct}}.  It also saves 
#'		this object as an .RData file so that it can be used in 
#'		future analyses.
#'		
#' @export
structure2conStruct <- function(infile,start.loci,missing.datum,outfile){
	structure.data <- utils::read.table(infile,header=FALSE,stringsAsFactors=FALSE)
	sample.names <- structure.data[,1]
	genos <- structure.data[,start.loci:ncol(structure.data)]
	rm(structure.data)
		if(ncol(genos) %% 2 != 0){
			stop("\nyou have mis-specified the genotype matrix\nplease check documentation\n\n")
		}
	n.loci <- ncol(genos)/2
	freqs <- get.freqs(genos,n.loci)
	missing.data <- get.missing.data(genos,n.loci,missing.datum)
	freqs[missing.data==2] <- NA
	row.names(freqs) <- sample.names
	outfile <- paste0(outfile,".RData")
	if(file.exists(outfile)){
		stop("\noutfile already exists\n\n")
	}
	save(freqs,file=outfile)
	return(freqs)
}

get.counted.allele <- function(genos){
	alleles <- unique(genos)
	alleles <- alleles[!alleles<0]
	counted <- sample(alleles,1)
	return(counted)
}

get.freqs <- function(genos,n.loci){
	freqs <- Reduce("cbind",
				lapply(1:n.loci,
							function(l){
								(genos[,seq(1,2*n.loci,by=2)[l]] == 1) + 
								(genos[,seq(2,2*n.loci,by=2)[l]] == 1)
							}))
	freqs <- freqs/2
	return(freqs)
}

get.missing.data <- function(genos,n.loci,missing.datum){
	missing.data <- Reduce("cbind",
						lapply(1:n.loci,
							function(l){
								(genos[,seq(1,2*n.loci,by=2)[l]] == missing.datum) + 
								(genos[,seq(2,2*n.loci,by=2)[l]] == missing.datum)
							}))
	return(missing.data)
}