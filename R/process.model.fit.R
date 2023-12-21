unstandardize.distances <- function(data.block){
	if(!is.null(data.block$sd.geoDist)){
		data.block$geoDist <- data.block$geoDist*data.block$sd.geoDist
	}
	return(data.block)
}

get.conStruct.results <- function(data.block,model.fit,n.chains){
	conStruct.results <- stats::setNames(
							lapply(1:n.chains,
								function(i){
									get.conStruct.chain.results(data.block,model.fit,i)
								}),
						  paste0("chain_",1:n.chains))
	return(conStruct.results)
}

get.MAP.iter <- function(model.fit,chain.no){
	lpd <- rstan::get_logposterior(model.fit)
	MAP.iter <- lapply(lpd,which.max)[[chain.no]]
	return(MAP.iter)
}

get.admix.props <- function(model.fit,chain.no,N,n.layers){
	# recover()
	admix.props <- array(1,dim=c(model.fit@sim$n_save[chain.no],N,n.layers))
	if(any(grepl("w",model.fit@model_pars))){
		for(k in 1:n.layers){
			admix.props[,,k] <- rstan::extract(model.fit,
											pars=unlist(lapply(1:N,function(j){sprintf("w[%s,%s]",j,k)})),
											permuted=FALSE,inc_warmup=TRUE)[,chain.no,]
		}
	}
	return(admix.props)
}

get.par.cov <- function(model.fit,chain.no,N){
	par.cov <- array(NA,dim=c(model.fit@sim$n_save[chain.no],N,N))
	for(i in 1:N){
		for(j in 1:N){
			my.par <- sprintf("parCov[%s,%s]",i,j)
			par.cov[,i,j] <- rstan::extract(model.fit,pars=my.par,inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
		}
	}
	return(par.cov)
}

get.nuggets <- function(model.fit,chain.no,N){
	nuggets <- rstan::extract(model.fit,pars="nugget",inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(nuggets)
}

get.gamma <- function(model.fit,chain.no){
	gamma <- rstan::extract(model.fit,pars="gamma",inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(gamma)
}

get.null.alpha.params <- function(n.iter){
	alpha.params <- list("alpha0" = rep(0,n.iter),
						 "alphaD" = rep(0,n.iter),
						 "alpha2" = rep(0,n.iter))
	return(alpha.params)	
}

get.alpha.params <- function(model.fit,data.block,chain.no,layer,n.layers){
	alpha.pars <- model.fit@model_pars[grepl("alpha",model.fit@model_pars)]
	if(length(alpha.pars) !=0 ){
		if(n.layers > 1){
			alpha.params <- stats::setNames(
									lapply(1:length(alpha.pars),
											function(i){
												rstan::extract(model.fit,
														pars=paste0(alpha.pars[i],"[",layer,"]"),
														inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
											}),alpha.pars)
		} else {
			alpha.params <- stats::setNames(
									lapply(1:length(alpha.pars),
											function(i){
												rstan::extract(model.fit,
														pars=alpha.pars[i],
														inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
											}),alpha.pars)		
		}
	} else {
		alpha.params <- get.null.alpha.params(model.fit@sim$n_save[chain.no])
	}
	if(!is.null(data.block$sd.geoDist)){
		alpha.params$alphaD <- alpha.params$alphaD/data.block$sd.geoDist
	}
	return(alpha.params)
}

get.null.phi <- function(n.iter){
	phi <- rep(0,n.iter)
	return(phi)
}

get.layer.phi <- function(model.fit,chain.no,layer){
	has.phi <- any(grepl("phi",model.fit@model_pars))
	if(has.phi){
		phi <- rstan::extract(model.fit,
						pars=paste0("phi","[",layer,"]"),
						inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	} else {
		phi <- get.null.phi(model.fit@sim$n_save[chain.no])
	}
	return(phi)
}

get.cov.function <- function(data.block){
	if(data.block$K == 1){
		if(data.block$spatial){
			cov.func <- function(layer.params,data.block){
				return(layer.params$alpha0 * 
						exp(-(layer.params$alphaD*data.block$geoDist)^layer.params$alpha2))
			}
		}
		if(!data.block$spatial){
			cov.func <- function(layer.params,data.block){
				return(matrix(0,nrow=data.block$N,ncol=data.block$N))
			}
		}
	} else {
		if(data.block$spatial){
			cov.func <- function(layer.params,data.block){
				return(layer.params$alpha0 *  
						exp(-(layer.params$alphaD*data.block$geoDist)^layer.params$alpha2) + 
							layer.params$phi)
			}
		}
		if(!data.block$spatial){
			cov.func <- function(layer.params,data.block){
				return(matrix(layer.params$phi,nrow=data.block$N,ncol=data.block$N))
			}
		}
	}
	return(cov.func)
}

get.layer.cov <- function(layer.params,data.block,n.iter){
	cov.function <- get.cov.function(data.block)
	layer.cov <- lapply(1:n.iter,
							function(i){
								cov.function(layer.params=
												lapply(layer.params,"[[",i),
												data.block)
							})
	return(layer.cov)
}

get.layer.params <- function(model.fit,data.block,chain.no,layer,n.layers,n.iter){
	layer.params <- list()
	layer.params <- get.alpha.params(model.fit,data.block,chain.no,layer,n.layers)
	layer.params[["phi"]] <- get.layer.phi(model.fit,chain.no,layer)
	layer.cov <- get.layer.cov(layer.params,data.block,n.iter)
	layer.params <- c(layer.params,list("layer.cov"=layer.cov))
	return(layer.params)
}

get.layer.params.list <- function(model.fit,data.block,chain.no,n.iter){
	layer.params <- stats::setNames(
								lapply(1:data.block$K,
											function(i){
												get.layer.params(model.fit,data.block,chain.no,i,data.block$K,n.iter)
											}),
								paste("layer",1:data.block$K,sep="_"))
	layer.params <- make.layer.params.S3(layer.params)
	return(layer.params)
}

make.layer.params.S3 <- function(layer.params){
	layer.params <- layer.params
	class(layer.params) <- "layer.params"
	return(layer.params)
}

#' An S3 print method for class layer.params
#' 
#' @param x an object of class \code{layer.params}
#' @param ... further options to be passed to \code{print}
#' @return prints a top-level summary of the layer.params, returns nothing
#' @method print layer.params
print.layer.params <- function(x,...){
	print(x=utils::str(x,max.level=1),...)
}

index.MAP <- function(param,MAP.iter){
	if(inherits(param,"numeric")){
		MAP.param <- param[MAP.iter]
	}
	if(inherits(param,"list")){
		MAP.param <- param[[MAP.iter]]
	}
	if(inherits(param,"array") & length(dim(param)) == 3){
		MAP.param <- param[MAP.iter,,]
		if(is.null(dim(MAP.param))){
			MAP.param <- matrix(MAP.param,nrow=length(MAP.param),ncol=1)
		}
	}
	if(inherits(param,"matrix") & length(dim(param)) == 2){
		MAP.param <- param[MAP.iter,]
	}
	if(inherits(param,"layer.params")){
		MAP.param <- index.MAP.layer.params.list(param,MAP.iter)
	}
	return(MAP.param)
}

index.MAP.layer.params <- function(layer.params,MAP.iter){
	MAP.layer.params <- lapply(layer.params,index.MAP,MAP.iter)
	return(MAP.layer.params)
}

index.MAP.layer.params.list <- function(layer.params.list,MAP.iter){
	MAP.layer.params.list <- lapply(layer.params.list,index.MAP.layer.params,MAP.iter)
	return(MAP.layer.params.list)
}

get.n.iter <- function(model.fit,chain.no){
	n.iter <- model.fit@sim$n_save[chain.no]
	return(n.iter)
}

make.conStruct.results.S3 <- function(conStruct.results){
	conStruct.results <- conStruct.results
	class(conStruct.results) <- "conStruct.results"
	return(conStruct.results)
}

#' An S3 print method for class conStruct.results
#' 
#' @param x an object of class \code{conStruct.results}
#' @param ... further options to be passed to \code{print}
#' @return prints a top-level summary of the conStruct.results, returns nothing
#' @method print conStruct.results
print.conStruct.results <- function(x,...){
	print(x=utils::str(x,max.level=1),...)
}

get.conStruct.chain.results <- function(data.block,model.fit,chain.no){
	n.iter <- get.n.iter(model.fit,chain.no)
	posterior <- list("n.iter" = model.fit@sim$n_save[chain.no],
					  "lpd" = rstan::get_logposterior(model.fit)[[chain.no]],
					  "nuggets" = get.nuggets(model.fit,chain.no,data.block$N),
					  "par.cov" = get.par.cov(model.fit,chain.no,data.block$N),
					  "gamma" = get.gamma(model.fit,chain.no),
					  "layer.params" = get.layer.params.list(model.fit,data.block,chain.no,n.iter),
					  "admix.proportions" = get.admix.props(model.fit,chain.no,data.block$N,data.block$K))
	MAP.iter <- get.MAP.iter(model.fit,chain.no)
	MAP <- lapply(posterior,function(X){index.MAP(X,MAP.iter)})
	names(MAP)[[1]]  <- "index.iter"
	MAP[["index.iter"]] <- MAP.iter
	conStruct.results <- list("posterior" = posterior,"MAP" = MAP)
	conStruct.results <- make.conStruct.results.S3(conStruct.results)
	return(conStruct.results)
}