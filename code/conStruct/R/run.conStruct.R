#' Run a conStruct analysis.
#'
#' \code{conStruct} runs a conStruct analysis of genetic data.
#'
#' This function initiates an analysis that uses  
#' geographic and genetic relationships between samples 
#' to estimate sample membership (admixture proportions) across 
#' a user-specified number of layers.
#'
#' @param spatial A logical indicating whether to perform a spatial analysis. 
#' 				  Default is \code{TRUE}. 
#' @param K An \code{integer} that indicates the number of layers to be 
#' 				  included in the analysis.
#' @param freqs A \code{matrix} of allele frequencies with one column per 
#'				locus and one row per sample.
#' 				Missing data should be indicated with \code{NA}.
#' @param geoDist A \code{matrix} of geographic distance between samples. 
#'					If \code{NULL}, user can only run the nonspatial model.
#' @param coords A \code{matrix} giving the longitude and latitude 
#'					(or X and Y coordinates) of the samples.
#' @param prefix A character \code{vector} giving the prefix to be attached 
#'					 to all output files.
#' @param n.chains An integer indicating the number of MCMC chains to be run 
#'					in the analysis. Default is 1.
#' @param n.iter An \code{integer} giving the number of iterations each MCMC 
#'				 chain is run. Default is 1e3.  If the number of iterations 
#'				 is greater than 500, the MCMC is thinned so that the number 
#'				 of retained iterations is 500 (before burn-in).
#' @param make.figs A \code{logical} value indicating whether to automatically 
#'					make figures once the analysis is complete. Default is 
#'					\code{TRUE}.
#' @param save.files A \code{logical} value indicating whether to automatically 
#'						save output and intermediate files once the analysis is
#'						 complete. Default is \code{TRUE}.
#'
#' @return This function returns a list with one entry for each chain run 
#'			(specified with \code{n.chains}). The entry for each chain is named 
#'			"chain_X" for the Xth chain.  The components of the entries for each 
#'			are detailed below: 
#'			\itemize{
#'				\item \code{posterior} gives parameter estimates over the posterior 
#'						distribution of the MCMC.
#'					\itemize{
#'						\item \code{n.iter} number of MCMC iterations retained for 
#'								analysis (half of the \code{n.iter} argument 
#'								specified in the function call).
#'						\item \code{lpd} vector of log posterior density over the retained 
#'								MCMC iterations.
#'						\item \code{nuggets} matrix of estimated nugget parameters with 
#'								one row per MCMC iteration and one column per sample.
#'						\item \code{par.cov} array of estimated parametric covariance matrices, 
#'								for which the first dimension is the number of MCMC iterations.
#'						\item \code{gamma} vector of estimated gamma parameter.
#'						\item \code{cluster.params} list summarizing estimates of layer-specific 
#'								parameters. There is one entry for each layer specified, and the 
#'								entry for the kth layer is named "Cluster_k".
#'							\itemize{
#'								\item \code{alpha0} vector of estimated alpha0 parameter in the 
#'										kth cluster.
#'								\item \code{alphaD} vector of estimated alphaD parameter in the 
#'										kth cluster.
#'								\item \code{alpha2} vector of estimated alpha2 parameter in the 
#'										kth cluster.
#'								\item \code{mu} vector of estimated mu parameter in the 
#'										kth cluster.
#'								\item \code{cluster.cov} vector of estimated cluster-specific 
#'										covariance parameter in the kth cluster.
#'							}
#'						\item \code{admix.proportions} array of estimated admixture proportions.
#'								The first dimension is the number of MCMC iterations, 
#'								the second is the number of samples, 
#' 								and the third is the number of layers.
#'					}
#'			\item \code{MAP} gives point estimates of the parameters listed in the \code{posterior}
#'								list described above. Values are indexed at the MCMC iteration 
#'								with the greatest posterior probability.
#'					\itemize{
#'						\item \code{index.iter} the iteration of the MCMC with the highest 
#'								posterior probability, which is used to index all parameters 
#'								included in the \code{MAP} list
#'						\item \code{lpd} the greatest value of the posterior probability
#'						\item \code{nuggets} point estimate of nugget parameters
#'						\item \code{par.cov} point estimate of parametric covariance
#'						\item \code{gamma} point estimate of gamma parameter
#'						\item \code{cluster.params} point estimates of all cluster-specific parameters 
#'						\item \code{admix.proportions} point estimates of admixture proportions.
#'					}
#'			}
#'
#' @details This function acts as a wrapper around a STAN model block determined 
#'			by the user-specified model (e.g., a spatial model with 3 layers, 
#'			or a nonspatial model with 5 layers).
#'			User-specified data are checked for appropriate format and consistent dimensions,
#'			then formatted into a \code{data.block},
#'			which is then passed to the STAN model block.
#'			Along with the \code{conStruct.results} output described above, 
#'			several objects are saved during the course of a \code{conStruct} call
#'			(if \code{save.files=TRUE}).
#'			These are the \code{data.block}, which contains all data passed to the STAN model block,
#'			\code{model.fit}, which is unprocessed results of the STAN run in \code{stanfit} format,
#'			and the \code{conStruct.results}, which are saved in the course of the function call
#'			in addition to being returned.
#'			If \code{make.figs=TRUE}, running \code{conStruct} will also generate many output figures, 
#'			which are detailed in the function \code{make.all.the.plots} in this package.
#'
#' @import rstan
#' @export
conStruct <- function(spatial=TRUE,K,freqs,geoDist=NULL,coords,prefix="",n.chains=1,n.iter=1e3,make.figs=TRUE,save.files=TRUE){
	call.check <- check.call(args <- as.list(environment()))
	freq.data <- process.freq.data(freqs)
	data.block <- make.data.block(K,freq.data,coords,spatial,geoDist,temp=NULL)
		if(save.files){
			save(data.block,file=paste0(prefix,"_data.block.Robj"))
		}
	stan.block <- make.stan.code.block(spatial,K)
	model.fit <- rstan::stan(model_code = stan.block,
							 refresh = min(n.iter/10,500),
							 data = data.block,
							 iter = n.iter,
							 chains = n.chains,
							 thin = ifelse(n.iter/500 > 1,n.iter/500,1),
							 save_warmup = FALSE)
	#save fit obj
		if(save.files){
			save(model.fit,file=paste(prefix,"model.fit.Robj",sep="_"))
		}
	conStruct.results <- get.conStruct.results(data.block,model.fit,n.chains)
		if(save.files){
			save(conStruct.results,file=paste(prefix,"conStruct.results.Robj",sep="_"))
		}
	if(make.figs){
		make.all.the.plots(conStruct.results,data.block,prefix,cluster.colors=NULL)
	}
	return(conStruct.results)
}

validate.data.list <- function(data.block){
	if(!"spatial" %in% names(data.block)){
		stop("\nUser must specify a \"spatial\" option\n\n")
	}
	if(!"N" %in% names(data.block)){
		stop("\nUser must specify a \"N\"\n\n")
	}
	if(!"K" %in% names(data.block)){
		stop("\nUser must specify a \"K\"\n\n")
	}
	if(!"L" %in% names(data.block)){
		stop("\nUser must specify a \"L\"\n\n")
	}
	if(!"obsCov" %in% names(data.block)){
		stop("\nUser must specify a \"obsCov\"\n\n")
	}
	return(invisible("list elements validated"))	
}

validate.n.samples <- function(data.block){
	n.samples <- data.block$N
	n.samples <- c(data.block$N,nrow(data.block$obsCov))
	if(!is.null(data.block$geoDist)){
		n.samples <- c(n.samples,nrow(data.block$geoDist))
	}
	if(length(unique(n.samples)) > 1){
		stop("\nthe number of samples is not consistent 
				across entries in the data.block\n\n")
	}
	return(invisible("n.samples validated"))
}

validate.model <- function(data.block){
	if(data.block$spatial){
		if(is.null(data.block$geoDist)){
			stop("\nyou have specified a spatial model,
				  but you have not specified a matrix 
				  of pairwise geographic distances\n\n")
		}
		if(any(data.block$geoDist < 0)){
			stop("\nyou have specified an invalid 
				  distance matrix that contains 
				  negative values\n\n")
		}
		if(any(is.na(data.block$geoDist))){
			stop("\nyou have specified an invalid 
				  distance matrix that contains 
				  non-numeric values\n\n")			
		}
	}
	return(invisible("model validated"))
}

make.data.block.S3 <- function(data.block){
	data.block <- data.block
	class(data.block) <- "data.block"
	return(data.block)
}

print.data.block <- function(data.block){
	print(utils::str(data.block,max.level=1))
}

validate.data.block <- function(data.block){
	message("\nchecking data.block\n")
		validate.data.list(data.block)
		validate.n.samples(data.block)
	message(sprintf("\treading %s samples",data.block$N))
	message(sprintf("\treading %s loci",data.block$L))
	message("\nchecking specified model\n")
		validate.model(data.block)
	if(data.block$spatial){
		message(sprintf("\nuser has specified a spatial model with %s cluster(s)\n",data.block$K))
	}
	if(!data.block$spatial){
		message(sprintf("\nuser has specified a purely discrete model with %s cluster(s)\n",data.block$K))
	}
	data.block <- make.data.block.S3(data.block)
	return(data.block)
}

make.stan.code.block <- function(spatial,n.clusters){
	stan.code.block.name <- "stan.block"
	if(n.clusters == 1){
		stan.code.block.name <- paste0("oneK.",stan.code.block.name)
	}
	if(n.clusters > 1){
		stan.code.block.name <- paste0("multiK.",stan.code.block.name)
	}
	if(spatial){
		stan.code.block.name <- paste0("space.",stan.code.block.name)
	}
	return(get(stan.code.block.name))
}

make.freq.data.list.S3 <- function(freq.data){
	freq.data <- freq.data
	class(freq.data) <- "freq.data"
	return(freq.data)
}

print.freq.data <- function(freq.data){
	print(utils::str(freq.data,max.level=1))
}

identify.invar.sites <- function(freqs){
	invar <- length(unique(freqs[which(!is.na(freqs))])) == 1
	return(invar)
}

drop.invars <- function(freqs){
	invars <- apply(freqs,2,identify.invar.sites)
	freqs <- freqs[,!invars]
	return(freqs)
}

identify.missing.sites <- function(freqs){
	n.samples <- length(freqs)
	missing <- FALSE
	if(length(which(is.na(freqs))) == n.samples){
		missing <- TRUE
	}
	return(missing)
}

drop.missing <- function(freqs){
	missing <- apply(freqs,2,identify.missing.sites)
	freqs <- freqs[,!missing]
	return(freqs)
}

calc.covariance <- function(freqs){
	x <- t(freqs)
	allelic.covariance <- stats::cov(x,use="pairwise.complete.obs") - 
									(1/2) * outer( colMeans(x,na.rm=TRUE), 1-colMeans(x,na.rm=TRUE), "*" ) -
									(1/2) * outer(1-colMeans(x,na.rm=TRUE), colMeans(x,na.rm=TRUE), "*") + 1/4
	diag(allelic.covariance) <- 1/4
	return(allelic.covariance)
}

process.freq.data <- function(freqs){
	freqs <- drop.invars(freqs)
	freqs <- drop.missing(freqs)
	n.loci <- ncol(freqs)
	obsCov <- calc.covariance(freqs)
	freq.data <- list("freqs" = freqs,
					  "obsCov" = obsCov,
					  "n.loci" = n.loci)
	freq.data <- make.freq.data.list.S3(freq.data)
	return(freq.data)
}

standardize.distances <- function(D){
	if(!is.null(D)){
		stdev.D <- stats::sd(D[upper.tri(D)])
		std.D <- D/stdev.D
	} else {
		std.D <- NULL
	}
	return(std.D)
}

set.temp <- function(temp=NULL){
	if(is.null(temp)){
		temp <- 1
	} else {
		temp <- temp
	}
	return(temp)
}

make.data.block <- function(K,freq.data,coords,spatial,geoDist=NULL,temp=NULL){
	data.block <- list("N" = nrow(coords),
					   "K" = K,
					   "spatial" = spatial,
					   "L" = freq.data$n.loci,
					   "coords" = coords,
					   "obsCov" = freq.data$obsCov,
					   "geoDist" = standardize.distances(geoDist),
					   "varMeanFreqs" = mean(0.5*colMeans(freq.data$freqs-0.5,na.rm=TRUE)^2 + 0.5*colMeans(1-freq.data$freqs-0.5,na.rm=TRUE)^2),
					   "temp" = set.temp(temp))
	data.block <- validate.data.block(data.block)
	return(data.block)
}

check.call <- function(args){
	if(args[["spatial"]] != TRUE & args[["spatial"]] != FALSE){
		stop("\nyou have specified an invalid value for the \"spatial\" argument \n")
	}
	if(length(args[["K"]]) > 1 | class(args[["K"]]) == "character"){
		stop("\nyou have specified an invalid value for the \"K\" argument \n")
	}
	if(class(args[["freqs"]]) != "matrix" | any(args[["freqs"]] > 1,na.rm=TRUE) | any(args[["freqs"]] < 0,na.rm=TRUE)){	
		stop("\nyou have specified an invalid value for the \"freqs\" argument \n")
	}
	if(class(args[["geoDist"]]) != "matrix" | length(unique(dim(args[["geoDist"]]))) > 1 | any(args[["D"]] < 0)){	
		stop("\nyou have specified an invalid value for the \"geoDist\" argument \n")
	}
	if(class(args[["coords"]]) != "matrix" | ncol(args[["coords"]]) > 2){	
		stop("\nyou have specified an invalid value for the \"coords\" argument \n")
	}
	return(invisible("args checked"))		
}