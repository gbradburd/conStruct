#'@export
conStruct <- function(spatial=TRUE,K,freqs,geoDist=NULL,temp=NULL,coords,prefix="",n.chains=1,n.iter=1e3,make.figs=TRUE,save.files=TRUE){
	call.check <- check.call(args <- as.list(environment()))
	#validate data block
	freq.data <- process.freq.data(freqs)
	data.block <- make.data.block(K,freq.data,coords,spatial,geoDist,temp)
		if(save.files){
			save(data.block,file=paste0(prefix,"_data.block.Robj"))
		}
	#validate model specification
	#make.stan.code.block
	stan.block <- make.stan.code.block(spatial,K)
		#write stan block to file
	#run model
	#put stan in tryCatch, email me
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
		make.all.the.plots(conStruct.results,n.chains,data.block,prefix,cluster.colors=NULL)
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
	print(str(data.block,max.level=1))
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
	print(str(freq.data,max.level=1))
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
	allelic.covariance <- cov(x,use="pairwise.complete.obs") - 
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
		stdev.D <- sd(D[upper.tri(D)])
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