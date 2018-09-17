#' Run a conStruct cross-validation analysis
#' 
#' \code{x.validation} runs a conStruct cross-validation analysis
#' 
#' This function initiates a cross-validation analysis that 
#' uses Monte Carlo cross-validation to determine the statistical 
#' support for models with different numbers of layers or 
#' with and without a spatial component.
#' 
#' @param train.prop A numeric value between 0 and 1 that gives 
#'		the proportions of the data to be used in the 
#'		training partition of the analysis. Default is 0.9.
#' @param n.reps An \code{integer} giving the number of cross-
#'		validation replicates to be run.
#' @param K A numeric \code{vector} giving the numbers of layers 
#'		to be tested in each cross-validation replicate.
#'		E.g., \code{K=1:7}.
#' @param freqs A \code{matrix} of allele frequencies with one column per 
#'		locus and one row per sample.
#' 		Missing data should be indicated with \code{NA}.
#' @param data.partitions A list with one element for each desired 
#'		cross-validation replicate. This argument can be specified 
#'		instead of the \code{freqs} argument if the user wants to 
#'		provide their own data partitions for model training and testing.
#'		See the model comparison vignette for details on what this 
#'		should look like.
#' @param geoDist A \code{matrix} of geographic distance between samples. 
#'		If \code{NULL}, user can only run the nonspatial model.
#' @param coords A \code{matrix} giving the longitude and latitude 
#'		(or X and Y coordinates) of the samples.
#' @param prefix A character \code{vector} giving the prefix to be attached 
#'		to all output files.
#' @param n.iter An \code{integer} giving the number of iterations each MCMC 
#'		chain is run. Default is 1e3.  If the number of iterations 
#'		is greater than 500, the MCMC is thinned so that the number 
#'		of retained iterations is 500 (before burn-in).
#' @param make.figs A \code{logical} value indicating whether to automatically 
#'		make figures during the course of the cross-validation analysis. 
#'		Default is \code{FALSE}.
#' @param save.files A \code{logical} value indicating whether to automatically 
#'		save output and intermediate files once the analysis is
#'		complete. Default is \code{FALSE}.
#' @param parallel A \code{logical} value indicating whether or not to run the 
#'		different cross-validation replicates in parallel. Default is \code{FALSE}.
#'		For more details on how to set up runs in parallel, see the model 
#'		comparison vignette.
#' @param n.nodes Number of nodes to run parallel analyses on. Default is 
#'		\code{NULL}. Ignored if \code{parallel} is \code{FALSE}. For more details 
#'		in how to set up runs in parallel, see the model comparison vignette. 
#'	
#' @return This function returns (and also saves as a .Robj) a \code{list} 
#'		containing the standardized results of the cross-validation analysis
#'		across replicates.  For each replicate, the function returns 
#' 		a list with the following elements:
#' 	\itemize{
#' 		\item \code{sp} - the mean of the standardized log likelihoods of the 
#'		"testing" data partition of that replicate for the spatial model for
#' 		each value of K specified in \code{K}.
#' 		\item \code{nsp} - the mean of the standardized log likelihoods of the 
#'		"testing" data partitions of that replicate for the nonspatial model for
#' 		each value of K specified in \code{K}.
#' }
#' In addition, this function saves two text files containing the standardized 
#'	cross-validation results for the spatial and nonspatial results 
#'	(prefix_sp_xval_results.txt and prefix_nsp_xval_results.txt, respectively).
#'	These values are written as matrices for user convenience; each column is 
#'	a cross-validation replicate, and each row gives the result for a value of 
#'	\code{K}.
#'
#'@export
x.validation <- function(train.prop = 0.9, n.reps, K, freqs = NULL, data.partitions = NULL, geoDist, coords, prefix, n.iter, make.figs = FALSE, save.files = FALSE,parallel=FALSE,n.nodes=NULL){
	call.check <- check.xval.call(args <- as.list(environment()))
	if(is.null(data.partitions)){
		data.partitions <- make.data.partitions(n.reps,freqs,train.prop)
	}
	check.data.partitions.arg(args <- as.list(environment()))
	save(data.partitions,file=paste0(prefix, ".xval.data.partitions.Robj"))
	prespecified <- parallel.prespecify.check(args <- as.list(environment()))
	`%d%` <- parallelizing(args <- as.list(environment()))
	i <- 1
    x.val <- foreach::foreach(i=1:n.reps) %d% {
        				x.validation.rep(rep.no = i, 
        								 K, 
        								 data.partition = data.partitions[[i]], 
        								 geoDist, 
        								 coords, 
        								 prefix, 
        								 n.iter, 
        								 make.figs, 
        								 save.files)
    			 }
    names(x.val) <- paste0("rep_", 1:n.reps)
    x.val <- lapply(x.val, standardize.xvals)
    save(x.val,file=paste0(prefix,".xval.results.Robj"))
	write.xvals(x.val,prefix)
	tmp <- end.parallelization(prespecified)
    return(x.val)
}

#' Match layers up across independent conStruct runs
#' 
#' \code{match.layers.x.runs} 
#' 
#' This function takes the results of two independent
#' \code{conStruct} analyses and compares them to identify 
#' which layers in a new analysis correspond most closely 
#' to the layers from an original analysis.
#' 
#' @param admix.mat1 A \code{matrix} of estimated admixture proportions
#'			from the original \code{conStruct} analysis, with one row 
#'			per sample and one column per layer. 
#' @param admix.mat2 A \code{matrix} of estimated admixture proportions
#'			from a second \code{conStruct} analysis, with one row per 
#'			sample and one column per layer, for which the 
#'			layer order is desired. Must have equal or greater number 
#'			of layers to \code{admix.mat1}.
#' @param admix.mat1.order An optional \code{vector} giving the  
#'			order in which the layers of \code{admix.mat1} are read.
#' 
#' @return This function returns a \code{vector} giving the ordering 
#' 			of the layers in \code{admix.mat2} that maximizes 
#' 			similarity between \code{admix.mat1} and re-ordered 
#'			\code{admix.mat2}.
#' 
#' @details This function compares admixture proportions in layers across 
#'			independent \code{conStruct} runs, and compares between them to 
#' 			identify the layers in \code{admix.mat2} that correspond most 
#'			closely to those in \code{admix.mat1}. It then returns a vector 
#' 			giving an ordering of \code{admix.mat2} that matches up the order
#' 			of the layers that correspond to each other.  This can be useful 
#' 			for:
#'			\enumerate{
#'				\item Dealing with "label switching" across independent runs 
#'					with the same number of layers; 
#'				\item Plotting results from independent runs with different 
#'					numbers of layers using consistent colors
#' 					(e.g., the "blue" layer shows up as blue even as 
#'					\code{K} increases); 
#' 				\item Examining results for multimodality (i.e., multiple 
#'					distinct solutions with qualitatively different patterns
#'					of membership across layers).
#' 			}
#' 			The \code{admix.mat1.order} argument can be useful when running 
#' 			this function to sync up plotting colors/order across the output 
#' 			of more than two \code{conStruct} runs.
#' 
#' @examples
#' \dontshow{
#'		admix.props1 <- matrix(c(0.09,0.00,0.50,0.51,0.10,0.05,0.02,0.01,0.80,0.00,0.22,0.74,0.92,0.20,0.47,0.00,0.78,0.30,0.33,0.45,0.00,0.00,0.64,0.90,0.00,0.00,0.00,0.01,0.02,0.00,0.00,0.09,0.00,0.55,0.00,0.00,0.00,0.09,0.02,0.00,0.00,0.01,0.00,0.20,0.00,0.06,0.05,0.08,0.04,0.01,0.00,0.06,0.17,0.14,0.03,0.00,0.00,0.18,0.08,0.00,1.00,1.00,0.99,0.98,0.98,1.00,0.74,0.98,0.43,1.00,0.91,1.00,0.41,0.47,0.90,0.95,0.96,0.99,0.00,1.00,0.72,0.20,0.00,0.77,0.52,1.00,0.15,0.53,0.53,0.53,1.00,1.00,0.18,0.02,1.00,0.00,0.00,0.00,0.00,0.02,0.00,0.17,0.02,0.01,0.00),ncol=3)
#'		admix.props2 <- matrix(c(0.36,0.35,0.42,0.38,0.35,0.35,0.36,0.35,0.48,0.36,0.39,0.39,0.40,0.36,0.36,0.35,0.40,0.46,0.45,0.38,0.34,0.35,0.47,0.40,0.35,1.00,1.00,0.99,0.99,0.98,1.00,0.84,0.99,0.63,1.00,0.32,0.35,0.24,0.24,0.33,0.34,0.33,0.35,0.15,0.32,0.32,0.10,0.30,0.33,0.27,0.36,0.13,0.26,0.27,0.22,0.36,0.35,0.14,0.11,0.35,0.00,0.00,0.00,0.01,0.01,0.00,0.07,0.00,0.18,0.00,0.32,0.30,0.34,0.38,0.31,0.30,0.31,0.30,0.36,0.32,0.30,0.51,0.30,0.31,0.37,0.30,0.47,0.29,0.28,0.40,0.30,0.31,0.39,0.49,0.30,0.00,0.00,0.00,0.00,0.01,0.00,0.09,0.01,0.19,0.00),ncol=3)
#' }
#' # compare the estimated admixture proportions from 
#'	# two different conStruct runs to determine which 
#'	# layers in one run correspond to those in the other
#' match.layers.x.runs(admix.props1,admix.props2)
#' 
#'@export
match.layers.x.runs <- function(admix.mat1,admix.mat2,admix.mat1.order=NULL){
	#recover()
	K1 <- ncol(admix.mat1)
		if(!is.null(admix.mat1.order)){
			admix.mat1 <- admix.mat1[,admix.mat1.order]
		}
	K2 <- ncol(admix.mat2)
	if(K1 > K2){
		stop("\nadmix.mat1 cannot have more layers than admix.mat2\n")
	}
	k.combn <- expand.grid(1:K1,1:K2)
	layer.sims <- unlist(lapply(1:nrow(k.combn),
						function(n){
							measure.frob.similarity(admix.mat1[,k.combn[n,1],drop=FALSE],
													admix.mat2[,k.combn[n,2],drop=FALSE],
													K=1)
							}))
	run2.order <- numeric(K2)
	while(length(which(run2.order == 0)) > (K2-K1)){
		tmp.max <- which.max(rank(layer.sims,na.last=FALSE))
		run2.match <- k.combn[tmp.max,2]
		run1.match <- k.combn[tmp.max,1]
		run2.order[run2.match] <- run1.match
		layer.sims[which(k.combn[,1]==run1.match)] <- NA
		layer.sims[which(k.combn[,2]==run2.match)] <- NA
	}
	if(K2 > K1){
		run2.order[which(run2.order==0)] <- (K1+1):K2
	}
	run2.order <- order(run2.order)
	return(run2.order)
}

#' Calculate layer contribution
#' 
#' \code{calculate.layer.contribution} 
#' 
#' This function takes the results of a \code{conStruct} 
#' analysis and calculates the relative contributions of 
#' each layer to total covariance.
#' 
#' @param conStruct.results The list output by a 
#'			\code{conStruct} run for a given MCMC chain. 
#' @param data.block A \code{data.block} list saved during a 
#'			\code{conStruct} run.
#' @param layer.order An optional \code{vector} giving the  
#'			order in which the layers of \code{conStruct.results} are 
#' 			read.
#' 
#' @return This function returns a \code{vector} giving the 
#' 			relative contributions of the layers 
#' 			in the analysis.
#' 
#' @details This function calculates the contribution of each layer to
#'			total covariance by multiplying the within-layer covariance 
#'			in a given layer by the admixture proportions samples draw 
#'			from that layer. The relative contribution of that layer 
#'			is this absolute contribution divided by the sum of those of 
#'			all other layers. 
#' 			A layer can have a large contribution if many samples draw 
#'			large amounts of admixture from it, or if it has a very large 
#'			within-layer covariance parameter (phi), or some combination 
#'			of the two. Layer contribution can be useful for evaluating 
#'			an appropriate level of model complexity for the data (e.g., 
#'			choosing a value of \code{K} or comparing the spatial and 
#'			nonspatial models).
#'@export
calculate.layer.contribution <- function(conStruct.results,data.block,layer.order=NULL){
	if(any(grepl("chain",names(conStruct.results)))){
		stop("user must specify conStruct results from a single chain\ni.e. from conStruct.results[[1]] rather than conStruct.results")
	}
	K <- ncol(conStruct.results$MAP$admix.proportions)
	apply.over <- 1:K
	if(!is.null(layer.order)){
		apply.over <- apply.over[layer.order]
	}
	raw.layer.scores <- lapply(apply.over,function(k){
							calculate.laycon.k(k,conStruct.results$MAP,data.block)
						})
	layer.contributions <- unlist(raw.layer.scores)/sum(unlist(raw.layer.scores))
	return(layer.contributions)
}

make.data.partitions <- function(n.reps,freqs,train.prop){
	data.partitions <- lapply(1:n.reps,
							function(i){
								xval.process.data(freqs,train.prop)								
							})
	names(data.partitions) <- paste0("rep",1:n.reps)
	return(data.partitions)
}

get.var.mean.freqs <- function(freqs){
	varMeanFreqs <- mean(0.5 * colMeans(freqs - 0.5, na.rm = TRUE)^2 + 
	            		 0.5 * colMeans(0.5 - freqs, na.rm = TRUE)^2)
	return(varMeanFreqs)
}

xval.process.data <- function(freqs,train.prop){
    freqs <- drop.invars(freqs)
    train.loci <- sample(1:ncol(freqs), ncol(freqs) * (train.prop))
    test.loci <- c(1:ncol(freqs))[!(1:ncol(freqs)) %in% train.loci]
    train.data <- calc.covariance(freqs[, train.loci])
    test.data <- calc.covariance(freqs[, test.loci])
    data.partition <- list("training" = list("data" = train.data,
    										 "n.loci" = length(train.loci),
    										 "varMeanFreqs" = get.var.mean.freqs(freqs[,train.loci])),
						   "testing" = list("data" = test.data,
    										 "n.loci" = length(test.loci),
    										 "varMeanFreqs" = get.var.mean.freqs(freqs[,test.loci])))
    return(data.partition)
}

xval.make.data.block <- function(K, data.partition, coords, spatial, geoDist = NULL){
	sd.dist.list <- standardize.distances(geoDist)
	data.block <- list(N = nrow(coords),
					   K = K,
					   spatial = spatial,
					   L = data.partition$n.loci,
					   coords = coords,
					   obsCov = data.partition$data,
					   geoDist = sd.dist.list$std.D,
					   sd.geoDist = sd.dist.list$stdev.D,
					   varMeanFreqs = data.partition$varMeanFreqs)
    data.block <- validate.data.block(data.block)
    return(data.block)
}

check.data.partitions.covmats <- function(args){
	data.list1 <- lapply(args[["data.partitions"]],function(x){x[[1]]$data})
	data.list2 <- lapply(args[["data.partitions"]],function(x){x[[2]]$data})
	if(!all(unlist(lapply(data.list1,function(x){isSymmetric(x)}))) |
	   !all(unlist(lapply(data.list2,function(x){isSymmetric(x)})))){
		stop("\nyou must specify symmetric matrices for the \"data\" elemtns of the data partitions list\n\n")
	}	
	if(any(unlist(lapply(data.list1,function(x){any(is.na(x))}))) |
	   any(unlist(lapply(data.list2,function(x){any(is.na(x))})))){
		stop("\nyou have specified an invalid data partition \"data\" element that contains non-numeric elements\n\n")
	}
	if(any(unlist(lapply(data.list1,function(x){any(eigen(x)$values <= 0)}))) |
	   any(unlist(lapply(data.list2,function(x){any(eigen(x)$values <= 0)})))){
		stop("\nyou have specified an invalid data partition \"data\" element is not positive definite\n\n")
	}
	return(invisible("data partitions cov matrices checked"))
}

check.data.partitions.arg <- function(args){
	if(length(args[["data.partitions"]]) != args[["n.reps"]]){
		stop("\nyou must specify 1 data partition for each cross-validation rep\n\n")
	}
	if(!all(unlist(lapply(args[["data.partitions"]],function(x){names(x)==c("training","testing")})))){
		stop("\neach element of the data partition list must contain an element named \"training\" and an element named \"testing\"\n\n")
	}
	if(!all(unlist(lapply(args[["data.partitions"]],function(x){names(x[[1]])==c("data","n.loci","varMeanFreqs")}))) | 
	   !all(unlist(lapply(args[["data.partitions"]],function(x){names(x[[1]])==c("data","n.loci","varMeanFreqs")})))){
		stop("\nwithin each element of the data partition named \"training\" or \"testing\", there must be three elements named \"data\",\"n.loci\", and \"varMeanFreqs\"\n\n")
	}
	L.list <- unlist(lapply(args[["data.partitions"]],function(x){c(x[[1]]$n.loci,x[[2]]$n.loci)}))
	vmf.list <- unlist(lapply(args[["data.partitions"]],function(x){c(x[[1]]$varMeanFreqs,x[[2]]$varMeanFreqs)}))
	if(any(is.na(L.list)) | any(L.list < 0) | any(!is.numeric(L.list))){
		stop("\nall values of \"n.loci\" must contain a numeric value greater than zero\n\n")			
	}
	if(any(is.na(vmf.list)) | any(vmf.list < 0) | any(!is.numeric(vmf.list))){
		stop("\nall values of \"varMeanFreqs\" must contain a numeric value greater than zero\n\n")			
	}
	check.data.partitions.covmats(args)
	return(invisible("data partitions arg checked"))
}

check.genetic.data.arg <- function(args){
	if(is.null(args[["data.partitions"]]) & is.null(args[["freqs"]])){
		stop("\nyou must specify a value for either the \"freqs\" argument or the \"data.partitions\" argument\n\n")
	}
	if(is.null(args[["data.partitions"]])){
		check.freqs.arg(args)
	}
	if(is.null(args[["freqs"]])){
		check.data.partitions.arg(args)
	}
	return(invisible("genetic data args checked"))
}

check.for.files <- function(args){
	if(file.exists(paste0(args[["prefix"]],"_sp_xval_results.txt")) |
	   file.exists(paste0(args[["prefix"]],"_nsp_xval_results.txt")) |
	   file.exists(paste0(args[["prefix"]],".xval.results.Robj"))){
		stop("\noutput files will be overwritten if you proceed with this analysis\n\n")
	}
	return(invisible("files checked for"))
}

check.parallel.args <- function(args){
	if(args[["parallel"]] & args[["n.nodes"]]==1){
		stop("\nyou have specified the \"parallel\" option with \"n.nodes\" set to 1.\n\n")
	}
	if(!args[["parallel"]] & args[["n.nodes"]] > 1){
		stop("\nyou have are running with \"parallel\" set to FALSE but with \"n.nodes\" greater than 1.\n\n")
	}
	if(!args[["parallel"]] & foreach::getDoParWorkers() > 1){
		stop("\nyou are running with more than 1 worker but you have set the \"parallel\" option to FALSE\n\n")
	}
	return(invisible("parallel args checked"))
}

check.xval.call <- function(args){
	check.for.files(args)
	check.genetic.data.arg(args)
	args$spatial <- TRUE
	check.geoDist.arg(args)
	check.coords.arg(args)
	check.parallel.args(args)
	return(invisible("args checked"))
}

xval.conStruct <- function (spatial = TRUE, K, data, geoDist = NULL, coords, prefix = "", n.chains = 1, n.iter = 1000, make.figs = TRUE, save.files = TRUE) {
    data.block <- xval.make.data.block(K, data, coords, spatial, geoDist)
    if (save.files) {
        save(data.block, file = paste0(prefix, "_data.block.Robj"))
    }
	stan.model <- pick.stan.model(spatial,K)
    model.fit <- rstan::sampling(object = stanmodels[[stan.model]], 
    							 refresh = min(n.iter/10,500), 
    							 data = data.block, 
    							 iter = n.iter, 
    							 chains = n.chains, 
        						 thin = ifelse(n.iter/500 > 1, n.iter/500, 1), 
        						 save_warmup = FALSE)
    conStruct.results <- get.conStruct.results(data.block,model.fit,n.chains)
	data.block <- unstandardize.distances(data.block)
    if (save.files) {
        save(data.block, file = paste0(prefix, "_data.block.Robj"))	
        save(conStruct.results, file = paste(prefix, "conStruct.results.Robj", sep = "_"))
        save(model.fit, file = paste(prefix, "model.fit.Robj", sep = "_"))
    }
    if (make.figs) {
        make.all.the.plots(conStruct.results, data.block, prefix,layer.colors = NULL)
    }
    return(conStruct.results)
}

x.validation.rep <- function(rep.no, K, data.partition, geoDist, coords, prefix, n.iter, make.figs = FALSE, save.files = FALSE) {
    training.runs.sp <- lapply(K, function(k) {
        xval.conStruct(spatial = TRUE, K = k, 
					   data = data.partition$training, 
					   geoDist = geoDist, coords = coords, 
					   prefix = paste0(prefix, "_sp_", "rep", rep.no, "K", k), 
					   n.iter = n.iter, make.figs = make.figs, save.files = save.files)
    })
    names(training.runs.sp) <- paste0("K", K)
    if (save.files) {
        save(training.runs.sp, file = paste0(prefix, "_rep", 
            rep.no, "_", "training.runs.sp.Robj"))
    }
    training.runs.nsp <- lapply(K, function(k) {
        xval.conStruct(spatial = FALSE, K = k, 
					   data = data.partition$training, 
					   geoDist = geoDist, coords = coords, 
					   prefix = paste0(prefix, "_nsp_", "rep", rep.no, "K", k), 
					   n.iter = n.iter, make.figs = make.figs, save.files = save.files)
    })
    names(training.runs.nsp) <- paste0("K", K)
    if (save.files) {
        save(training.runs.nsp, file = paste0(prefix, "_rep", 
            rep.no, "_", "training.runs.nsp.Robj"))
    }
    test.lnl.sp <- lapply(training.runs.sp, function(x) {
        fit.to.test(data.partition$testing, x[[1]])
    })
    	names(test.lnl.sp) <- K
    test.lnl.nsp <- lapply(training.runs.nsp, function(x) {
        fit.to.test(data.partition$testing, x[[1]])
    })
    	names(test.lnl.nsp) <- K
    test.lnl <- list(sp = test.lnl.sp, nsp = test.lnl.nsp)
    if (save.files) {
        save(test.lnl, file = paste0(prefix, "rep", rep.no, "_test.lnl.Robj"))
    }
    return(test.lnl)
}

standardize.xvals <- function(x.val){
	mean.lnls <- lapply(x.val,function(x){
					lapply(x,function(k){mean(unlist(k))})
				 })
	xval.max <- max(unlist(mean.lnls))
	mean.std.lnls <- lapply(mean.lnls,function(s){
						lapply(s,function(k){
							k - xval.max
						})})
	return(mean.std.lnls)
}

get.xval.CIs <- function(x.vals.std,K){
	#recover()
	sp.means <- lapply(
				   		lapply(1:K,
							function(k){
								unlist(lapply(
									lapply(x.vals.std,"[[","sp"),
						"[[",k))}),
					function(x){
						mean(x)})
	sp.std.errs <- lapply(
				   		lapply(1:K,
							function(k){
								unlist(lapply(
									lapply(x.vals.std,"[[","sp"),
						"[[",k))}),
					function(x){
						stats::sd(x)/sqrt(length(x))})
	sp.CIs <- lapply(1:K,function(k){
					sp.means[[k]] + c(-1.96*sp.std.errs[[k]],
									   1.96*sp.std.errs[[k]])})
	nsp.means <- lapply(
				   		lapply(1:K,
							function(k){
								unlist(lapply(
									lapply(x.vals.std,"[[","nsp"),
						"[[",k))}),
					function(x){
						mean(x)})
	nsp.std.errs <- lapply(
				   		lapply(1:K,
							function(k){
								unlist(lapply(
									lapply(x.vals.std,"[[","nsp"),
						"[[",k))}),
					function(x){
						stats::sd(x)/sqrt(length(x))})
	nsp.CIs <- lapply(1:K,function(k){
					nsp.means[[k]] + c(-1.96*nsp.std.errs[[k]],
									    1.96*nsp.std.errs[[k]])})
	return(list("sp.means" = unlist(sp.means),
				"sp.std.errs" = unlist(sp.std.errs),
				"sp.CIs" = sp.CIs,
				"nsp.means" = unlist(nsp.means),
				"nsp.std.errs" = unlist(nsp.std.errs),
				"nsp.CIs" = nsp.CIs))
}

write.xvals <- function(xvals,prefix){
	sp <- data.frame(lapply(lapply(xvals,"[[","sp"),unlist))
		names(sp) <- paste0("sp_",names(sp))
		sp <- round(sp,digits=4)
	nsp <- data.frame(lapply(lapply(xvals,"[[","nsp"),unlist))
		names(nsp) <- paste0("nsp_",names(nsp))
		nsp <- round(nsp,digits=4)
	utils::write.table(sp,file=paste0(prefix,"_sp_xval_results.txt"),row.names=FALSE,quote=FALSE)
	utils::write.table(nsp,file=paste0(prefix,"_nsp_xval_results.txt"),row.names=FALSE,quote=FALSE)
}

parallel.prespecify.check <- function(args){
	prespecified <- FALSE
	if(args[["parallel"]] & foreach::getDoParRegistered()){
		prespecified <- TRUE
	}
	return(prespecified)
}

end.parallelization <- function(prespecified){
	if(!prespecified){
		doParallel::stopImplicitCluster()
		message("\nParallel workers terminated\n\n")
	}
	return(invisible("if not prespecified, parallelization ended"))
}

parallelizing <- function(args){
	if(args[["parallel"]]){
		if(!foreach::getDoParRegistered()){
			if(is.null(args[["n.nodes"]])){
				n.nodes <- parallel::detectCores()-1
			} else {
				n.nodes <- args[["n.nodes"]]
			}
			cl <- parallel::makeCluster(n.nodes)
			doParallel::registerDoParallel(cl)
			message("\nRegistered doParallel with ",n.nodes," workers\n")
		} else {
			message("\nUsing ",foreach::getDoParName()," with ", foreach::getDoParWorkers(), " workers")
		}
		d <- foreach::`%dopar%`
	} else {
		message("\nRunning sequentially with a single worker\n")
		d <- foreach::`%do%`
	}
	return(d)
}

calculate.qij <- function(layer.params,data.block,i,j){
	q_ij <- 2 * layer.params$alpha0 * 
				exp(-(layer.params$alphaD * 
						data.block$geoDist[i,j])^layer.params$alpha2) + 
						2 * layer.params$phi + 0.5
	return(q_ij)
}

calculate.weighted.qij <- function(k,MAP,data.block,i,j){
	#recover()
	w.qij <- MAP$admix.proportions[i,k] * MAP$admix.proportions[j,k] * 
				calculate.qij(MAP$layer.params[[k]],data.block,i,j)
	return(w.qij)
}

calculate.laycon.k <- function(k,MAP,data.block){
	comps <- cbind(utils::combn(1:data.block$N,2),sapply(1:data.block$N,rep,2))
	lay.con <- lapply(1:ncol(comps),function(n){
					calculate.weighted.qij(k,MAP,data.block,comps[1,n],comps[2,n])
				})
	return(mean(unlist(lay.con)))
}

post.process.par.cov <- function(conStruct.results,samples){
	pp.cov.list <- lapply(samples,
							function(i){
								list("inv" = chol2inv(chol(conStruct.results$posterior$par.cov[i,,])),
									 "log.det" = determinant(conStruct.results$posterior$par.cov[i,,])$modulus[[1]])
							})
	return(pp.cov.list)
}


log.likelihood <- function(obsCov,inv.par.cov,log.det,n.loci){
	lnL <- -0.5 * (sum( inv.par.cov * obsCov) + n.loci * log.det)
	return(lnL)
}

calc.lnl.x.MCMC <- function(cov.chunk,pp.par.cov,chunk.size){
	lnl.x.mcmc <- lapply(pp.par.cov,
						function(x){
							log.likelihood(cov.chunk,x$inv,x$log.det,n.loci=chunk.size)
					})
	return(unlist(lnl.x.mcmc))
}

fit.to.test <- function(test.data,conStruct.results){
	pp.par.cov <- post.process.par.cov(conStruct.results,
										samples = 1:conStruct.results$posterior$n.iter)
	test.lnl <- lapply(pp.par.cov,
						function(x){
							log.likelihood(test.data$data,x$inv,x$log.det,test.data$n.loci)
				})
	return(test.lnl)
}

frob.mn <- function(M){
	frob.mn <- sqrt(sum(M^2))
	return(frob.mn)
}

measure.frob.similarity <- function(m1,m2,K){
	W <- matrix(1/K,nrow(m1),ncol(m1))
	frob.sim <- 1 - frob.mn(m1-m2)/sqrt(frob.mn(m1-W) * frob.mn(m2-W))
	return(frob.sim)
}