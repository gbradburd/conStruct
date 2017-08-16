get.conStruct.results <- function(data.block,model.fit,n.chains){
	conStruct.results <- setNames(
						lapply(1:n.chains,
							function(i){
								get.conStruct.chain.results(data.block,model.fit,i)
							}),
					  paste0("chain_",1:n.chains))
	return(conStruct.results)
}

get.MAP.iter <- function(model.fit,chain.no){
	lpd <- get_logposterior(model.fit)
	MAP.iter <- lapply(lpd,which.max)[[chain.no]]
	return(MAP.iter)
}

get.admix.props <- function(model.fit,chain.no,N,n.clusters){
	# recover()
	admix.props <- array(1,dim=c(model.fit@sim$n_save[chain.no],N,n.clusters))
	if(any(grepl("w",model.fit@model_pars))){
		for(k in 1:n.clusters){
			admix.props[,,k] <- extract(model.fit,
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
			par.cov[,i,j] <- extract(model.fit,pars=my.par,inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
		}
	}
	return(par.cov)
}

get.nuggets <- function(model.fit,chain.no,N){
	nuggets <- extract(model.fit,pars="nugget",inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(nuggets)
}

get.gamma <- function(model.fit,chain.no){
	gamma <- extract(model.fit,pars="gamma",inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(gamma)
}

get.null.alpha.params <- function(n.iter){
	alpha.params <- list("alpha0" = rep(0,n.iter),
						 "alphaD" = rep(0,n.iter),
						 "alpha2" = rep(0,n.iter))
	return(alpha.params)	
}

get.alpha.params <- function(model.fit,chain.no,cluster,n.clusters){
	alpha.pars <- model.fit@model_pars[grepl("alpha",model.fit@model_pars)]
	if(length(alpha.pars) !=0 ){
		if(n.clusters > 1){
			alpha.params <- setNames(
								lapply(1:length(alpha.pars),
										function(i){
											extract(model.fit,
													pars=paste0(alpha.pars[i],"[",cluster,"]"),
													inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
										}),alpha.pars)
		} else {
			alpha.params <- setNames(
								lapply(1:length(alpha.pars),
										function(i){
											extract(model.fit,
													pars=alpha.pars[i],
													inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
										}),alpha.pars)		
		}
	} else {
		alpha.params <- get.null.alpha.params(model.fit@sim$n_save[chain.no])
	}
	return(alpha.params)
}

get.null.mu <- function(n.iter){
	mu <- rep(0,n.iter)
	return(mu)
}

get.cluster.mu <- function(model.fit,chain.no,cluster){
	has.mu <- any(grepl("mu",model.fit@model_pars))
	if(has.mu){
		mu <- extract(model.fit,
						pars=paste0("mu","[",cluster,"]"),
						inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	} else {
		mu <- get.null.mu(model.fit@sim$n_save[chain.no])
	}
	return(mu)
}

get.cov.function <- function(data.block){
	if(data.block$K == 1){
		if(data.block$spatial){
			cov.func <- function(cluster.params,data.block){
				return(cluster.params$alpha0 * 
						exp(-(cluster.params$alphaD*data.block$geoDist)^cluster.params$alpha2))
			}
		}
		if(!data.block$spatial){
			cov.func <- function(cluster.params,data.block){
				return(matrix(0,nrow=data.block$N,ncol=data.block$N))
			}
		}
	} else {
		if(data.block$spatial){
			cov.func <- function(cluster.params,data.block){
				return(cluster.params$alpha0 *  
						exp(-(cluster.params$alphaD*data.block$geoDist)^cluster.params$alpha2) + 
							cluster.params$mu)
			}
		}
		if(!data.block$spatial){
			cov.func <- function(cluster.params,data.block){
				return(matrix(cluster.params$mu,nrow=data.block$N,ncol=data.block$N))
			}
		}
	}
	return(cov.func)
}

get.cluster.cov <- function(cluster.params,data.block,n.iter){
	cov.function <- get.cov.function(data.block)
	cluster.cov <- lapply(1:n.iter,
							function(i){
								cov.function(cluster.params=
												lapply(cluster.params,"[[",i),
												data.block)
							})
	return(cluster.cov)
}

get.cluster.params <- function(model.fit,data.block,chain.no,cluster,n.clusters,n.iter){
	cluster.params <- list()
	cluster.params <- get.alpha.params(model.fit,chain.no,cluster,n.clusters)
	cluster.params[["mu"]] <- get.cluster.mu(model.fit,chain.no,cluster)
	cluster.cov <- get.cluster.cov(cluster.params,data.block,n.iter)
	cluster.params <- c(cluster.params,list("cluster.cov"=cluster.cov))
	return(cluster.params)
}

get.cluster.params.list <- function(model.fit,data.block,chain.no,n.iter){
	cluster.params <- setNames(
						lapply(1:data.block$K,
									function(i){
										get.cluster.params(model.fit,data.block,chain.no,i,data.block$K,n.iter)
									}),
						paste("Cluster",1:data.block$K,sep="_"))
	cluster.params <- make.cluster.params.S3(cluster.params)
	return(cluster.params)
}

make.cluster.params.S3 <- function(cluster.params){
	cluster.params <- cluster.params
	class(cluster.params) <- "cluster.params"
	return(cluster.params)
}

print.cluster.params <- function(cluster.params){
	print(str(cluster.params,max.level=1))
}

index.MAP <- function(param,MAP.iter){
	if(class(param) == "numeric"){
		MAP.param <- param[MAP.iter]
	}
	if(class(param) == "list"){
		MAP.param <- param[[MAP.iter]]
	}
	if(class(param) == "array"){
		MAP.param <- param[MAP.iter,,]
		if(is.null(dim(MAP.param))){
			MAP.param <- matrix(MAP.param,nrow=length(MAP.param),ncol=1)
		}
	}
	if(class(param) == "matrix"){
		MAP.param <- param[MAP.iter,]
	}
	if(class(param) == "cluster.params"){
		MAP.param <- index.MAP.cluster.params.list(param,MAP.iter)
	}
	return(MAP.param)
}

index.MAP.cluster.params <- function(cluster.params,MAP.iter){
	MAP.cluster.params <- lapply(cluster.params,index.MAP,MAP.iter)
	return(MAP.cluster.params)
}

index.MAP.cluster.params.list <- function(cluster.params.list,MAP.iter){
	MAP.cluster.params.list <- lapply(cluster.params.list,index.MAP.cluster.params,MAP.iter)
	return(MAP.cluster.params.list)
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

print.conStruct.results <- function(conStruct.results){
	print(str(conStruct.results,max.level=1))
}

get.conStruct.chain.results <- function(data.block,model.fit,chain.no){
	n.iter <- get.n.iter(model.fit,chain.no)
	posterior <- list("n.iter" = model.fit@sim$n_save[chain.no],
					  "lpd" = get_logposterior(model.fit)[[chain.no]],
					  "nuggets" = get.nuggets(model.fit,chain.no,data.block$N),
					  "par.cov" = get.par.cov(model.fit,chain.no,data.block$N),
					  "gamma" = get.gamma(model.fit,chain.no),
					  "cluster.params" = get.cluster.params.list(model.fit,data.block,chain.no,n.iter),
					  "admix.proportions" = get.admix.props(model.fit,chain.no,data.block$N,data.block$K))
	MAP.iter <- get.MAP.iter(model.fit,chain.no)
	MAP <- lapply(posterior,function(X){index.MAP(X,MAP.iter)})
	conStruct.results <- list("posterior" = posterior,"MAP" = MAP)
	conStruct.results <- make.conStruct.results.S3(conStruct.results)
	return(conStruct.results)
}
