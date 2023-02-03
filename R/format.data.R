#' Convert a dataset from STRUCTURE to conStruct format
#'
#' \code{structure2conStruct} converts a STRUCTURE dataset 
#' to conStruct format
#' 
#' This function takes a population genetics dataset in 
#' STRUCTURE format and converts it to conStruct format. 
#' The STRUCTURE file can have one row per individual 
#' and two columns per locus, or one column and two rows 
#' per individual. It can only contain bi-allelic SNPs.
#' Missing data is acceptable, but must be indicated with 
#' a single value throughout the dataset.
#' 
#' @param infile The name and path of the file in STRUCTURE format 
#' 			to be converted to \code{conStruct} format. 
#' @param onerowperind Indicates whether the file format has 
#'		one row per individual (\code{TRUE}) or two rows per 
#'		individual (\code{FALSE}).
#' @param start.loci The index of the first column in the dataset 
#'			that contains genotype data.
#' @param start.samples The index of the first row in the dataset 
#'			that contains genotype data (e.g., after any headers). 
#'			Default value is 1.
#' @param missing.datum The character or value used to denote 
#' 			missing data in the STRUCTURE dataset (often 0 or -9).
#' @param outfile The name and path of the file containing the 
#'			\code{conStruct} formatted dataset to be generated 
#' 			by this function.
#'
#' @details This function takes a STRUCTURE format data file and 
#'		converts it to a \code{conStruct} format data file.
#'		This function can only be applied to diploid organisms.
#'		The STRUCTURE data file must be a plain text file. 
#'		If there is extraneous text or column headers before the data 
#'		starts, those extra lines should be deleted by hand or 
#'		taken into account via the \code{start.samples} argument.
#'		
#' 		The STRUCTURE dataset can either be in the ONEROWPERIND=1 
#' 		file format, with one row per individual and two columns 
#' 		per locus, or the ONEROWPERIND=0 format, with two rows and 
#'		one column per individual. The first column of the STRUCTURE 
#' 		dataset should be individual names. There may be any number 
#' 		of other columns that contain non-genotype information before 
#'		the first column that contains genotype data, but there can 
#' 		be no extraneous columns at the end of the dataset, after the 
#' 		genotype data.
#'		
#'		The genotype data must be bi-allelic 
#'		single nucleotide polymorphisms (SNPs). Applying this function 
#'		to datasets with more than two alleles per locus may result in 
#'		cryptic failure. For more details, see the \code{format-data} 
#'		vignette.
#'	
#'	@return This function returns an allele frequency data matrix 
#'		that can be used as the \code{freqs} argument in a conStruct 
#'		analysis run using \code{\link{conStruct}}.  It also saves 
#'		this object as an .RData file so that it can be used in 
#'		future analyses.
#'		
#' @export
structure2conStruct <- function(infile,onerowperind,start.loci,start.samples=1,missing.datum,outfile){
	outfile <- paste0(outfile,".RData")
	if(file.exists(outfile)){
		stop("\noutfile already exists\n\n")
	}
	structure.data <- utils::read.table(infile,header=FALSE,skip=start.samples-1,stringsAsFactors=FALSE)
	sample.names <- get.sample.names(structure.data,onerowperind)
	genos <- structure.data[,start.loci:ncol(structure.data)]
	rm(structure.data)
	if(onerowperind & ncol(genos) %% 2 != 0){
		stop("\nyou have mis-specified the genotype matrix\nplease check documentation\n\n")
	}
	if(!onerowperind & nrow(genos) %% 2 != 0){
		stop("\nyou have mis-specified the genotype matrix\nplease check documentation\n\n")	
	}

	freqs <- get.freqs(genos,onerowperind,missing.datum)
	row.names(freqs) <- sample.names
	save(freqs,file=outfile)
	return(freqs)
}

get.sample.names <- function(structure.data,onerowperind){
	sample.names <- structure.data[,1]
	if(!onerowperind){
		sample.names <- sample.names[seq(1,length(sample.names),by=2)]
	}
	return(sample.names)
}

get.counted.allele <- function(genos,missing.datum){
	alleles <- unique(genos)
	if(all(alleles==missing.datum)){
		stop("\nyour dataset contains loci with all data missing. please remove and re-try.\n\n")
	}
	alleles <- alleles[!alleles==missing.datum]
	counted <- sample(alleles,1)
	return(counted)
}

get.freqs <- function(genos,onerowperind,missing.datum){
	n.loci <- ifelse(onerowperind,ncol(genos)/2,ncol(genos))
	if(onerowperind){
		freqs <- get.freqs.onerowperind(genos,n.loci,missing.datum)
	} else {
		freqs <- get.freqs.tworowperind(genos,n.loci,missing.datum)
	}
	colnames(freqs) <- NULL
	return(freqs)
}

get.freqs.onerowperind <- function(genos,n.loci,missing.datum){
	if(any(genos > 1)){
		counted.alleles <- apply(genos,2,get.counted.allele,missing.datum)
	} else {
		counted.alleles <- rep(1,n.loci)
	}
	freqs <- Reduce("cbind",
				lapply(1:n.loci,
							function(l){
								(genos[,seq(1,2*n.loci,by=2)[l]] == counted.alleles[l]) + 
								(genos[,seq(2,2*n.loci,by=2)[l]] == counted.alleles[l])
							}))
	freqs <- freqs/2
	missing.data <- Reduce("cbind",
						lapply(1:n.loci,
							function(l){
								(genos[,seq(1,2*n.loci,by=2)[l]] == missing.datum) + 
								(genos[,seq(2,2*n.loci,by=2)[l]] == missing.datum)
							}))
	freqs[missing.data==2] <- NA
	return(freqs)
}

get.freqs.tworowperind <- function(genos,n.loci,missing.datum){
	if(any(genos > 1)){
		counted.alleles <- apply(genos,2,get.counted.allele,missing.datum)
	} else {
		counted.alleles <- rep(1,n.loci)
	}
	freqs <- Reduce("cbind",
				lapply(1:n.loci,
							function(l){
								(genos[seq(1,nrow(genos),by=2),l] == counted.alleles[l]) + 
								(genos[seq(2,nrow(genos),by=2),l] == counted.alleles[l])
							}))
	freqs <- freqs/2
	missing.data <- Reduce("cbind",
						lapply(1:n.loci,
							function(l){
								(genos[seq(1,nrow(genos),by=2),l] == missing.datum) + 
								(genos[seq(2,nrow(genos),by=2),l] == missing.datum)
							}))
	freqs[missing.data==2] <- NA
	return(freqs)
}

