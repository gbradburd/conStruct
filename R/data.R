#' Example dataset used in a \code{conStruct} analysis
#' 
#' A simulated dataset containing the allele frequency 
#' and sampling coordinate data necessary to run a 
#' \code{conStruct} analysis.
#' 
#' @format A list with two elements:
#' \describe{
#'		\item{allele.frequencies}{a matrix with one row for each of 
#'			the 36 samples and one column for each of 10,000 loci, 
#'			giving the frequency of the counted allele at each locus 
#'			in each sample}
#'		\item{coords}{a matrix with one row for each of the 36 samples, 
#'			in the same order as that of the allele frequency matrix, 
#'			and two columns, the first giving the x-coordinate 
#'			(or longitude), the second giving the y-coordinate (or latitude)}
#' }
#
"conStruct.data"