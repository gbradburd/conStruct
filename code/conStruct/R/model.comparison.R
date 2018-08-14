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
#'
#' @return This function returns a \code{list} containing the 
#'		standardized results of the cross-validation analysis
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
#'@export
x.validation <- function(train.prop = 0.9, n.reps, K, freqs, geoDist, coords, prefix, n.iter, make.figs = FALSE, save.files = FALSE){
    x.val <- lapply(1:n.reps, 
    				function(i) {
        				x.validation.rep(rep.no = i, 
        								 train.prop, 
        								 K, 
        								 freqs, 
        								 geoDist, 
        								 coords, 
        								 prefix, 
        								 n.iter, 
        								 make.figs, 
        								 save.files)        				
    		 })
    names(x.val) <- paste0("rep_", 1:n.reps)
    x.val <- lapply(x.val, standardize.xvals)
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
#'			of the two. layer contribution can be useful for evaluating 
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


xval.process.data <- function(freqs,train.prop,rep.no,prefix){
    freqs <- drop.invars(freqs)
    train.loci <- sample(1:ncol(freqs), ncol(freqs) * (train.prop))
    test.loci <- c(1:ncol(freqs))[!(1:ncol(freqs)) %in% train.loci]
    train.data <- freqs[, train.loci]
    test.data <- freqs[, test.loci]
    xval.freq.data <- list("training" = train.data,
    					   "testing" = test.data)
    save(xval.freq.data, file = paste0(prefix, "_rep", rep.no, "_xval.dataset.Robj"))
    return(xval.freq.data)
}

xval.conStruct <- function (spatial = TRUE, K, freqs, geoDist = NULL, coords, prefix = "", n.chains = 1, n.iter = 1000, make.figs = TRUE, save.files = TRUE) {
    call.check <- check.call(args <- as.list(environment()))
    freq.data <- process.freq.data(freqs)
    data.block <- make.data.block(K, freq.data, coords, spatial, geoDist, temp = NULL)
    if (save.files) {
        save(data.block, file = paste0(prefix, "_data.block.Robj"))
    }
    model.fit <- rstan::sampling(object = pick.stan.model(spatial, K),
    							 refresh = min(n.iter/10,500), 
    							 data = data.block, 
    							 iter = n.iter, 
    							 chains = n.chains, 
        						 thin = ifelse(n.iter/500 > 1, n.iter/500, 1), 
        						 save_warmup = FALSE)
    if (save.files) {
        save(model.fit, file = paste(prefix, "model.fit.Robj", 
            sep = "_"))
    }
    conStruct.results <- get.conStruct.results(data.block,model.fit,n.chains)
    if (save.files) {
        save(conStruct.results, file = paste(prefix, "conStruct.results.Robj", 
            sep = "_"))
    }
    if (make.figs) {
        make.all.the.plots(conStruct.results, data.block, prefix,layer.colors = NULL)
    }
    return(conStruct.results)
}


x.validation.rep <- function(rep.no, train.prop, K, freqs, geoDist, coords, prefix, n.iter, make.figs = FALSE, save.files = FALSE) {
	xval.freq.data <- xval.process.data(freqs,train.prop,rep.no,prefix)
    training.runs.sp <- lapply(K, function(k) {
        xval.conStruct(spatial = TRUE, K = k, 
					   freqs = xval.freq.data$training, 
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
					   freqs = xval.freq.data$training, 
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
        fit.to.test(xval.freq.data$testing, x[[1]])
    })
    	names(test.lnl.sp) <- K
    test.lnl.nsp <- lapply(training.runs.nsp, function(x) {
        fit.to.test(xval.freq.data$testing, x[[1]])
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
	freq.data <- process.freq.data(test.data)
	pp.par.cov <- post.process.par.cov(conStruct.results,
										samples = 1:conStruct.results$posterior$n.iter)
	test.lnl <- lapply(pp.par.cov,
						function(x){
							log.likelihood(freq.data$obsCov,x$inv,x$log.det,freq.data$n.loci)
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
