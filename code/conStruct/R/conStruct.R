#'@export
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

#'@export
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

#'@export
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

#'@export
make.data.block.S3 <- function(data.block){
	data.block <- data.block
	class(data.block) <- "data.block"
	return(data.block)
}

#'@export
print.data.block <- function(data.block){
	print(str(data.block,max.level=1))
}

#'@export
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

#'@export
make.stan.code.block <- function(spatial,n.clusters){
	stan.code.block.name <- "stan.block"
	if(n.clusters == 1){
		stan.code.block.name<- paste0("oneK.",stan.code.block.name)
	}
	if(n.clusters > 1){
		stan.code.block.name<- paste0("multiK.",stan.code.block.name)
	}
	if(spatial){
		stan.code.block.name <- paste0("space.",stan.code.block.name)
	}
	return(get(stan.code.block.name))
}

#'@export
make.freq.data.list.S3 <- function(freq.data){
	freq.data <- freq.data
	class(freq.data) <- "freq.data"
	return(freq.data)
}

#'@export
print.freq.data <- function(freq.data){
	print(str(freq.data,max.level=1))
}

#'@export
haploidize <- function(freqs){
	pop.hap.var <- var(rbinom(length(freqs),1,freqs),na.rm=TRUE)
	return(pop.hap.var)
}

#'@export
unfold.freqs <- function(freqs){
	N <- nrow(freqs)
	L <- ncol(freqs)
	unfold <- runif(L,0,1) > 0.5
	freqs <- sapply(1:L,function(l){
				ifelse(rep(unfold[l],N),
						freqs[,l,drop=FALSE],
						1-freqs[,l,drop=FALSE])
			 })
	return(freqs)
}

#'@export
identify.invar.sites <- function(freqs){
	invar <- length(unique(freqs[which(!is.na(freqs))])) == 1
	return(invar)
}

#'@export
drop.invars <- function(freqs){
	invars <- apply(freqs,2,identify.invar.sites)
	freqs <- freqs[,!invars]
	return(freqs)
}

#'@export
identify.missing.sites <- function(freqs){
	n.samples <- length(freqs)
	missing <- FALSE
	if(length(which(is.na(freqs))) == n.samples){
		missing <- TRUE
	}
	return(missing)
}

#'@export
drop.missing <- function(freqs){
	missing <- apply(freqs,2,identify.missing.sites)
	freqs <- freqs[,!missing]
	return(freqs)
}

#'@export
process.freq.data <- function(freqs){
	freqs <- unfold.freqs(freqs)
	freqs <- drop.invars(freqs)
	freqs <- drop.missing(freqs)
	n.loci <- ncol(freqs)
	obsCov <- cov(t(freqs),use="pairwise.complete.obs")
	diag(obsCov) <- apply(freqs,1,haploidize)
	freq.data <- list("freqs" = freqs,
					  "obsCov" = obsCov,
					  "n.loci" = n.loci)
	freq.data <- make.freq.data.list.S3(freq.data)
	return(freq.data)
}


#'@export
standardize.distances <- function(D){
	if(!is.null(D)){
		stdev.D <- sd(D[upper.tri(D)])
		std.D <- D/stdev.D
	} else {
		std.D <- NULL
	}
	return(std.D)
}

#'@export
set.temp <- function(temp=NULL){
	if(is.null(temp)){
		temp <- 1
	} else {
		temp <- temp
	}
	return(temp)
}

#'@export
make.data.block <- function(K,freq.data,coords,spatial,geoDist=NULL,temp=NULL){
	data.block <- list("N" = nrow(coords),
					   "K" = K,
					   "spatial" = spatial,
					   "L" = freq.data$n.loci,
					   "coords" = coords,
					   "obsCov" = freq.data$obsCov,
					   "geoDist" = standardize.distances(geoDist),
					   "varMeanFreqs" = var(apply(freq.data$freqs,2,function(x){mean(x,na.rm=TRUE)})),
					   "temp" = set.temp(temp))
	data.block <- validate.data.block(data.block)
	return(data.block)
}

#'@export
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

#'@export
geoStructure <- function(spatial=TRUE,K,freqs,geoDist=NULL,temp=NULL,coords,prefix="",n.chains=1,n.iter=1e3,burnin=0){
	#MAKE IT SO BURNIN CAN'T BE MORE THAN N.ITER
	#GENERALLY DEAL W/ BURNIN
	call.check <- check.call(args <- as.list(environment()))
	#validate data block
	freq.data <- process.freq.data(freqs)
		save(freq.data,file=paste0(prefix,"_freq.data.Robj"))
	data.block <- make.data.block(K,freq.data,coords,spatial,geoDist,temp)
		save(data.block,file=paste0(prefix,"_data.block.Robj"))
	#validate model specification
	#make.stan.code.block
	stan.block <- make.stan.code.block(spatial,K)
		#write stan block to file
	#run model
	#put stan in tryCatch, email me
	require(rstan)
	model.fit <- stan(model_code = stan.block,
						refresh = min(n.iter/10,500),
						data = data.block,
						iter = n.iter,
						chains = n.chains,
						thin = ifelse(n.iter/500 > 1,n.iter/500,1))
	#save fit obj
	save(model.fit,file=paste(prefix,"model.fit.Robj",sep="_"))
	geoStr.results <- get.geoStructure.results(data.block,model.fit,n.chains)
	save(geoStr.results,file=paste(prefix,"geoStr.results.Robj",sep="_"))
	make.all.the.plots(geoStr.results,n.chains,data.block,freq.data,prefix,burnin,cluster.colors=NULL)
	return(geoStr.results)
}

#'@export
get.geoStructure.results <- function(data.block,model.fit,n.chains,burnin=0.2){
	geoStr.results <- setNames(
						lapply(1:n.chains,
							function(i){
								get.geoStructure.chain.results(data.block,model.fit,i,burnin)
							}),
					  paste0("chain_",1:n.chains))
	return(geoStr.results)
}

#'@export
get.MAP.iter <- function(model.fit,chain.no,burnin=0.5){
	burnin <- max(0.5,burnin)
	logpost <- get_logposterior(model.fit)
	MAP.iter <- lapply(logpost,
						function(x){
							length(x)*burnin + 
							which.max(x[(length(x)*burnin+1):length(x)])
						})[[chain.no]]
	return(MAP.iter)
}

#'@export
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

#'@export
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

#'@export
get.nuggets <- function(model.fit,chain.no,N){
	nuggets <- extract(model.fit,pars="nugget",inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(nuggets)
}

#'@export
get.gamma <- function(model.fit,chain.no){
	gamma <- extract(model.fit,pars="gamma",inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(gamma)
}

#'@export
get.alpha.params <- function(model.fit,chain.no,cluster,n.clusters){
	alpha.pars <- model.fit@model_pars[grepl("alpha",model.fit@model_pars)]
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
	return(alpha.params)
}

#'@export
get.cluster.mu <- function(model.fit,chain.no,cluster){
	mu <- extract(model.fit,
					pars=paste0("mu","[",cluster,"]"),
					inc_warmup=TRUE,permuted=FALSE)[,chain.no,]
	return(mu)
}

#'@export
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

#'@export
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

#'@export
get.cluster.params <- function(model.fit,data.block,chain.no,cluster,n.clusters,n.iter){
	cluster.params <- list()
	if(data.block$spatial){
		cluster.params <- get.alpha.params(model.fit,chain.no,cluster,n.clusters)
	}
	if(n.clusters > 1){
		cluster.params <- c(cluster.params,
							list("mu" = get.cluster.mu(model.fit,chain.no,cluster)))
	}
	cluster.cov <- get.cluster.cov(cluster.params,data.block,n.iter)
	cluster.params <- c(cluster.params,list("cluster.cov"=cluster.cov))
	return(cluster.params)
}

#'@export
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

#'@export
make.cluster.params.S3 <- function(cluster.params){
	cluster.params <- cluster.params
	class(cluster.params) <- "cluster.params"
	return(cluster.params)
}

#'@export
print.cluster.params <- function(cluster.params){
	print(str(cluster.params,max.level=1))
}

#'@export
index.MAP <- function(param,MAP.iter){
	if(class(param) == "numeric"){
		MAP.param <- param[MAP.iter]
	}
	if(class(param) == "list"){
		MAP.param <- param[[MAP.iter]]
	}
	if(class(param) == "array"){
		MAP.param <- param[MAP.iter,,]
	}
	if(class(param) == "matrix"){
		MAP.param <- param[MAP.iter,]
	}
	if(class(param) == "cluster.params"){
		MAP.param <- index.MAP.cluster.params.list(param,MAP.iter)
	}
	return(MAP.param)
}

#'@export
index.MAP.cluster.params <- function(cluster.params,MAP.iter){
	MAP.cluster.params <- lapply(cluster.params,index.MAP,MAP.iter)
	return(MAP.cluster.params)
}

#'@export
index.MAP.cluster.params.list <- function(cluster.params.list,MAP.iter){
	MAP.cluster.params.list <- lapply(cluster.params.list,index.MAP.cluster.params,MAP.iter)
	return(MAP.cluster.params.list)
}

#'@export
get.n.iter <- function(model.fit,chain.no){
	n.iter <- model.fit@sim$n_save[chain.no]
	return(n.iter)
}

#'@export
make.geoStructure.results.S3 <- function(geoStructure.results){
	geoStructure.results <- geoStructure.results
	class(geoStructure.results) <- "geoStructure.results"
	return(geoStructure.results)
}

#'@export
print.geoStructure.results <- function(geoStructure.results){
	print(str(geoStructure.results,max.level=1))
}

#'@export
get.geoStructure.chain.results <- function(data.block,model.fit,chain.no,burnin=0.2){
	n.iter <- get.n.iter(model.fit,chain.no)
	posterior <- list("n.iter" = model.fit@sim$n_save[chain.no],
					  "lpd" = get_logposterior(model.fit)[[chain.no]],
					  "nuggets" = get.nuggets(model.fit,chain.no,data.block$N),
					  "par.cov" = get.par.cov(model.fit,chain.no,data.block$N),
					  "gamma" = get.gamma(model.fit,chain.no))
	if(data.block$spatial | data.block$K > 1){
		posterior[["cluster.params"]] <- get.cluster.params.list(model.fit,data.block,chain.no,n.iter)
	}
	if(data.block$K > 1){
		posterior[["admix.proportions"]] <- get.admix.props(model.fit,chain.no,data.block$N,data.block$K)
	}
	MAP.iter <- get.MAP.iter(model.fit,chain.no,burnin)
	MAP <- lapply(posterior,function(X){index.MAP(X,MAP.iter)})
	geoStructure.results <- list("posterior" = posterior,"MAP" = MAP)
	geoStructure.results <- make.geoStructure.results.S3(geoStructure.results)
	return(geoStructure.results)
}

#'@export
plot.lpd <- function(geoStr.results,burnin=0){
	z <- (burnin+1):geoStr.results$posterior$n.iter
	x.lab <- ifelse(burnin==0,
						"MCMC iterations",
						paste0("MCMC iterations\n(post burnin of ",burnin,")"))
	plot(geoStr.results$posterior$lpd[z],
			xlab=x.lab,ylab="posterior probability",
			main="Posterior probability",type='l')
	return(invisible(0))
}

#'@export
plot.nuggets <- function(geoStr.results,burnin){
	z <- (burnin+1):geoStr.results$posterior$n.iter
	x.lab <- ifelse(burnin==0,
						"MCMC iterations",
						paste0("MCMC iterations\n(post burnin of ",burnin,")"))
	matplot(geoStr.results$post$nuggets[z,],type='l',
				main="sample nuggets",
				ylab="nugget value",
				xlab=x.lab)
	return(invisible("nuggets"))
}

plot.gamma <- function(geoStr.results,burnin){
	z <- (burnin+1):geoStr.results$posterior$n.iter
	x.lab <- ifelse(burnin==0,
						"MCMC iterations",
						paste0("MCMC iterations\n(post burnin of ",burnin,")"))
	plot(geoStr.results$posterior$gamma[z],
			xlab=x.lab,ylab="gamma",
			main="Gamma",type='l')
	return(invisible(0))
}

#'@export
get.ylim <- function(cluster.params,n.clusters,param,z){
	y.lim <- range(unlist(
				lapply(
					lapply(1:n.clusters,
						function(i){
							cluster.params[[i]][[param]]
						}),
					function(x){
						range(x[z])
					})))
	y.lim <- y.lim + c(-0.15*diff(y.lim),0.15*diff(y.lim))
	return(y.lim)
}

#'@export
plot.cluster.param <- function(cluster.param,clst.col,z){
	points(cluster.param[z],type='l',col=clst.col)
	return(invisible(0))
}

#'@export
plot.cluster.cov.params <- function(data.block,geoStr.results,burnin,cluster.colors){
	n.clusters <- data.block$K
	z <- (burnin+1):geoStr.results$posterior$n.iter
	x.lab <- ifelse(burnin==0,
						"MCMC iterations",
						paste0("MCMC iterations\n(post burnin of ",burnin,")"))
	params <- names(geoStr.results$posterior$cluster.params$Cluster_1)[!names(geoStr.results$posterior$cluster.params$Cluster_1)=="cluster.cov"]
	param.ranges <- lapply(params,function(x){get.ylim(geoStr.results$posterior$cluster.params,n.clusters,x,z)})
	if(length(params) > 0){
		for(i in 1:length(params)){
			plot(0,type='n',main=params[i],
				xlab=x.lab,ylab="parameter value",
				xlim=c(1,length(z)),
				ylim=param.ranges[[i]])
			lapply(1:n.clusters,function(j){plot.cluster.param(geoStr.results$posterior$cluster.params[[j]][[params[i]]],cluster.colors[j],z)})
			legend(x="topright",col= cluster.colors[1:n.clusters],lty=1,legend=paste0("Cluster_",1:n.clusters))
		}
	}
	return(invisible(0))
}

#'@export
plot.admix.props <- function(data.block,geoStr.results,cluster.colors,burnin){
	n.clusters <- data.block$K
	par(mfrow=c(n.clusters,1),mar=c(3,3,2,2))
		for(i in 1:n.clusters){
			matplot(geoStr.results$posterior$admix.proportions[(burnin+1):geoStr.results$posterior$n.iter,,i],type='l',ylim=c(0,1),
					main=paste0("Cluster ",i),ylab="admixture proportion",col=cluster.colors[i])
		}
	return(invisible(0))
}

#'@export
get.par.cov.CI <- function(data.block,geoStr.results,burnin){
	combns <- gtools::combinations(n=data.block$N,r=2,v=1:data.block$N,repeats.allowed=TRUE)
	iters <- (burnin+1):geoStr.results$posterior$n.iter
	CIs <- lapply(1:nrow(combns),
				function(i){
					quantile(geoStr.results$posterior$par.cov[iters,combns[i,1],combns[i,2]],c(0.025,0.975))
				})
	return(CIs)
}
#'@export
plot.model.fit.CIs <- function(data.block,freq.data,geoStr.results,burnin){
	cov.range <- range(c(freq.data$obsCov,
						geoStr.results$posterior$par.cov[
							(burnin+1):geoStr.results$posterior$n.iter, , ]))
	plot(data.block$geoDist,freq.data$obsCov,
    	xlab = "geographic distance", 
        ylab = "covariance",
        main="Cov/geoDist",
        ylim = cov.range, type = "n")
	combns <- gtools::combinations(n=data.block$N,r=2,v=1:data.block$N,repeats.allowed=TRUE)
	CIs <- get.par.cov.CI(data.block,geoStr.results,burnin)
	lapply(1:nrow(combns),
			function(i){
				segments(x0 = data.block$geoDist[combns[i,1],combns[i,2]],
						 y0 = CIs[[i]][1],
						 x1 = data.block$geoDist[combns[i,1],combns[i,2]],
						 y1 = CIs[[i]][2],
						 col = adjustcolor(1,0.1),
						 lwd=1.5)
			})
	points(data.block$geoDist,freq.data$obsCov,col=2,pch=20,cex=0.8)
	legend(x="topright",legend=c("observed","95% CI"),pch=c(19,NA),lty=c(NA,1),col=c(2,"gray"))
	return(invisible("plotted"))
}

#'@export
plot.model.fit <- function(data.block,freq.data,geoStr.results,burnin){
	z <- seq((burnin+1),geoStr.results$posterior$n.iter,length.out=10)
	index.mat <- upper.tri(data.block$geoDist, diag = TRUE)
	cov.range <- range(c(freq.data$obsCov,geoStr.results$posterior$par.cov[z, , ]))
    plot(data.block$geoDist,freq.data$obsCov,
    	xlab = "geographic distance", 
        ylab = "covariance",
        main="Cov/geoDist",
        ylim = cov.range, type = "n")
    lapply(z, function(i) {
        points(data.block$geoDist[index.mat], geoStr.results$posterior$par.cov[i,,][index.mat],
        	pch = 20, col = adjustcolor(1, 0.1))
    		})
    points(data.block$geoDist[index.mat], freq.data$obsCov[index.mat], 
        xlab = "geographic distance", ylab = "covariance", ylim = cov.range, 
        col=2,pch = 19)
	legend(x="topright",legend=c("observed","parametric"),pch=19,col=c(2,1))
	return(invisible("plotted"))
}

#'@export
plot.cluster.covs <- function(data.block,geoStr.results,cluster.colors,burnin){
	order.mat <- order(data.block$geoDist)
	z <- seq((burnin+1),geoStr.results$posterior$n.iter,length.out=10)
	y.range <- range(c(
				unlist(lapply(1:data.block$K,
							function(k){
								geoStr.results$posterior$cluster.params[[k]]$cluster.cov[z]})),
				data.block$obsCov))
	plot(data.block$geoDist[upper.tri(data.block$obsCov,diag=TRUE)],
		 data.block$obsCov[upper.tri(data.block$obsCov,diag=TRUE)],
			xlim=range(data.block$geoDist),ylim=y.range,
			xlab = "geographic distance",
			ylab = "cluster-specific covariances",
			main = "cluster covariances")
		lapply(1:data.block$K,function(k){
			lapply(geoStr.results$posterior$cluster.params[[k]]$cluster.cov[z],function(x){
				lines(data.block$geoDist[order.mat],
					  x[order.mat],col=adjustcolor(cluster.colors[k],0.3),pch=20)
			})
		})
	legend(x="topright",col= cluster.colors[1:data.block$K],lty=1,legend=paste0("Cluster_",1:data.block$K))
	return(invisible("cluster covs"))	
}

#'@export
plot.cluster.covariances <- function(data.block,geoStr.results,cluster.colors){
	ind.mat <- upper.tri(data.block$geoDist,diag=TRUE)
	y.range <- range(unlist(lapply(seq_along(1:data.block$K),
										function(i){geoStr.results$MAP$cluster.params[[i]]$cluster.cov})))
	plot(data.block$geoDist[ind.mat],data.block$geoDist[ind.mat],
			ylab="sample covariance",xlab="geographic distance",
			ylim=y.range + diff(range(y.range))/10 * c(-1,1),type='n')
	for(i in 1:data.block$K){
		points(data.block$geoDist[ind.mat],
				geoStr.results$MAP$cluster.params[[i]]$cluster.cov[ind.mat],
				col=cluster.colors[i],pch=20,cex=0.6)
	}
	legend(x="topright",lwd=2,lty=1,col=cluster.colors[1:data.block$K],
			legend=paste0("Cluster_",1:data.block$K),cex=0.7)
	return(invisible("plotted"))
}

#'@export
structure.polygon <- function(plotting.admix.props,i,j,use.colors){
	polygon(x = c(j-1,j,j,j-1),
			y = c(plotting.admix.props[i,j],
					plotting.admix.props[i,j],
					plotting.admix.props[i+1,j],
					plotting.admix.props[i+1,j]),
			col=use.colors[i])
	return(invisible(j))
}

#'@export
make.structure.polygon.layer <- function(plotting.admix.props,i,use.colors,sample.order){
		lapply(1:ncol(plotting.admix.props),function(j){
			structure.polygon(plotting.admix.props[,sample.order],i,j,use.colors)
		})
	return(invisible(i))
}

#'@export
make.structure.plot <- function(data.block,geoStr.results,mar=c(2,4,2,2),sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=NULL){
	par(mar=mar)
	if(is.null(cluster.order)){
		cluster.order <- seq(1:data.block$K)
	}
	if(is.null(sample.order)){
		sample.order <- seq(1:data.block$N)
	}
	if(!is.null(sort.by)){
		sample.order <- order(geoStr.results$MAP$admix.proportions[,sort.by])
	}
	if(is.null(cluster.colors)){
		cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	}
	if(data.block$K==1){
		geoStr.results$MAP$admix.proportions <- matrix(geoStr.results$MAP$admix.proportions,nrow=data.block$N,ncol=1)
	}
	use.colors <- cluster.colors[1:data.block$K][cluster.order]
	plot(0,xlim=c(0,data.block$N),ylim=c(0,1),type='n',ylab="admixture",xlab="",xaxt='n')
	plotting.admix.props <- apply(cbind(0,geoStr.results$MAP$admix.proportions[,cluster.order]),1,cumsum)
	lapply(1:data.block$K,function(i){
		make.structure.polygon.layer(plotting.admix.props,i,use.colors,sample.order)
	})
	if(!is.null(sample.names)){
		axis(side=1,at=seq(1:data.block$N)-0.5,labels=sample.names[sample.order],cex.axis=0.5,las=2)
	}
	return(invisible("plotted"))
}

#'@export
make.admix.pie.plot <- function(data.block,geoStr.results,cluster.colors,stat,radii=2.7,add=FALSE,title=NULL,x.lim=NULL,y.lim=NULL){
	if(is.null(data.block$coords)){
		message("\nuser has not specified sampling coordinates in the data block\n")
	} else {
		require(caroline)
		cluster.names <- paste0("cluster_",1:data.block$K)
		sample.names <- paste0("sample_",1:data.block$N)
		color.tab <- nv(c(cluster.colors[1:data.block$K]),cluster.names)
		if(stat == "MAP"){
			admix.props <- geoStr.results$MAP$admix.proportions
		} else if(stat == "mean"){
			admix.props <- apply(geoStr.results$posterior$admix.proportions,c(2,3),mean)
		} else if(stat == "median"){
			admix.props <- apply(geoStr.results$posterior$admix.proportions,c(2,3),median)		
		}
		pie.list <- lapply(1:data.block$N,function(i){nv(admix.props[i,],cluster.names)})
		names(pie.list) <- sample.names
		if(add){
			par(new=TRUE)
		} else {
			par(mar=c(2,2,2,2))
		}
		if(is.null(title)){
			title <- "Admixture proportion map"
		}
		if(is.null(x.lim)){
			x.lim <- c(min(data.block$coords[,1]) - 1, max(data.block$coords[,1]) + 1)
		}
		if(is.null(y.lim)){
			y.lim <- c(min(data.block$coords[,2]) - 1, max(data.block$coords[,2]) + 1)
		}
		pies(pie.list,x0=data.block$coords[,1],y0=data.block$coords[,2],
					color.table=color.tab,border="black",radii=radii,
					xlab="",ylab="",main=title,lty=1,density=NULL,
					xlim = x.lim, ylim = y.lim)
		box(lwd=2)
	}
	return(invisible(0))
}

#'@export
get.cluster.order <- function(K,admix.props,ref.admix.props){
	K.combn <- expand.grid(1:K,1:K)
	mean.props <- lapply(1:K,function(i){
						apply(admix.props[,,i],2,
							function(x){mean(x)})
					})	
	admix.prop.cors <- unlist(
						lapply(1:nrow(K.combn),function(i){
							cor(mean.props[[K.combn[i,1]]],
							ref.admix.props[,K.combn[i,2]])}))
	matchups <- matrix(c(1:K,rep(NA,K)),nrow=K,ncol=2)
		colnames(matchups) <- c("sample","reference")
	while(any(is.na(matchups))){
		winner <- as.numeric(K.combn[which.max(admix.prop.cors),])
		matchups[winner[1],2] <- winner[2]
		admix.prop.cors[which.max(admix.prop.cors)] <- NA
	}
	return(matchups)
}

#'@export
get.n.cluster.cov.params <- function(geoStr.results){
	n.params <- length(names(geoStr.results$posterior$cluster.params$Cluster_1)[
					!names(geoStr.results$posterior$cluster.params$Cluster_1)=="cluster.cov"])
	return(n.params)
}

#'@export
make.all.chain.plots <- function(geoStr.results,chain.no,data.block,freq.data=NULL,prefix,burnin=0,cluster.colors,...){
	pdf(file=paste0(prefix,"_trace.plots.chain_",chain.no,".pdf"),...)
		plot.lpd(geoStr.results,burnin)
		plot.nuggets(geoStr.results,burnin)
		plot.gamma(geoStr.results,burnin)
		plot.cluster.cov.params(data.block,geoStr.results,burnin,cluster.colors)
		if(data.block$K > 1){
			plot.admix.props(data.block,geoStr.results,cluster.colors,burnin)
		}
	dev.off()
	if(!is.null(freq.data)){
		pdf(file=paste0(prefix,"_model.fit.chain_",chain.no,".pdf"),...)
			plot.model.fit(data.block,freq.data,geoStr.results,burnin)
		dev.off()
		pdf(file=paste0(prefix,"_model.fit.CIs.chain_",chain.no,".pdf"),...)
			plot.model.fit.CIs(data.block,freq.data,geoStr.results,burnin)
		dev.off()
	}
	if(data.block$spatial | data.block$K > 1){
		pdf(file=paste0(prefix,"_cluster.cov.curves.chain_",chain.no,".pdf"),width=5,height=5)
			plot.cluster.covs(data.block,geoStr.results,cluster.colors,burnin)
		dev.off()
	}
	if(data.block$K > 1){
		pdf(file=paste0(prefix,"_pie.map.chain_",chain.no,".pdf"),width=6,height=6)	
			make.admix.pie.plot(data.block,geoStr.results,cluster.colors,stat="MAP",radii=2.7,add=FALSE,title=NULL,x.lim=NULL,y.lim=NULL)
		dev.off()
		pdf(file=paste0(prefix,"_structure.plot.chain_",chain.no,".pdf"),width=10,height=5)
			make.structure.plot(data.block,geoStr.results,mar=c(2,4,2,2),sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=cluster.colors)
		dev.off()
	}
	return(invisible("made chain plots!"))
}

#'@export
make.all.the.plots <- function(geoStr.results,n.chains,data.block,freq.data=NULL,prefix,burnin=0,cluster.colors=NULL,...){
	if(is.null(cluster.colors)){
		cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	}
	lapply(1:n.chains,function(i){
		make.all.chain.plots(geoStr.results[[i]],chain.no=i,data.block,freq.data,prefix,burnin,cluster.colors,...)
	})
	return(invisible("made chain plots!"))
}

#'@export
post.process.par.cov <- function(geoStr.results,samples){
	pp.cov.list <- lapply(samples,
							function(i){
								list("inv" = chol2inv(chol(geoStr.results$posterior$par.cov[i,,])),
									 "log.det" = determinant(geoStr.results$posterior$par.cov[i,,])$modulus[[1]])
							})
	return(pp.cov.list)
}

#'@export
log.likelihood <- function(obsCov,inv.par.cov,log.det,n.loci){
	lnL <- -0.5 * (sum( inv.par.cov * obsCov) + n.loci * log.det)
	return(lnL)
}

#'@export
determine.log.shift <- function(chunk.lnls){
	if(diff(range(unlist(chunk.lnls))) > 700){
		message("the difference between the min and max lnLs may be inducing underflow")
	}
	return(max(unlist(chunk.lnls)))
}

#'@export
shift.chunk.lnls <- function(chunk.lnls,A,n.iter){
	shift.chunk.lnls <- log(sum(exp(chunk.lnls-A))) + A - log(n.iter)
	return(shift.chunk.lnls)
}

#'@export
calc.lnl.x.MCMC <- function(cov.chunk,pp.par.cov,chunk.size){
	#recover()
	lnl.x.mcmc <- lapply(pp.par.cov,
						function(x){
							log.likelihood(cov.chunk,x$inv,x$log.det,n.loci=chunk.size)
					})
	return(unlist(lnl.x.mcmc))
}

#'@export
data.chunkification <- function(n.loci,chunk.size){
	chunk.intervals <- seq(1,n.loci,by=chunk.size)
	n.chunks <- length(chunk.intervals)	
	chunk.indices <- cbind(chunk.intervals[1:(n.chunks-1)],chunk.intervals[2:n.chunks])
	return(chunk.indices)
}

#'@export
chunk.freq.data <- function(freqs,data.block,chunk.size){
	n.loci <- ncol(freqs)
	chunk.indices <- data.chunkification(n.loci,chunk.size)
	chunks <- lapply(1:nrow(chunk.indices),
					  function(i){
					  	cov(t(freqs[,chunk.indices[i,1]:	
					  				 chunk.indices[i,2],drop=FALSE]),
					  		use="pairwise.complete.obs")
					  })
	return(chunks)
}

#'@export
calculate.lpd <- function(chunk.lnls,n.iter){
	#recover()
	#subtract max lnL from log likelihood to avoid overflow
	A <- determine.log.shift(chunk.lnls)
	chunk.lnls.shifted <- lapply(chunk.lnls,function(x){shift.chunk.lnls(x,A,n.iter)})
	lpd <- sum(unlist(chunk.lnls.shifted))
	return(lpd)
}

#'@export
calculate.pwaic <- function(chunk.lnls){
	pwaic <- lapply(chunk.lnls,var)
	return(sum(unlist(pwaic)))
}

#'@export
calculate.waic <- function(freqs,data.block,geoStr.results,chunk.size,samples=NULL){
	#recover()
	cat("breaking data into locus-by-locus covariances...\n\n")
	chunks <- chunk.freq.data(freqs,data.block,chunk.size)
	if(is.null(samples)){
		samples <- 1:geoStr.results$posterior$n.iter
	}
	n.iter <- length(samples)
	#invert the posterior distn of parametric cov matrices
	cat("inverting posterior distribution of parametric covariance matrices...\n\n")
		pp.par.cov <- post.process.par.cov(geoStr.results,samples)
	#calc likelioods of nth data chunk across all sampled MCMC iterations
	cat("calculating likelihood of each site across posterior distribution of parameters...\n\n")
		chunk.lnls <- lapply(chunks,function(x){calc.lnl.x.MCMC(x,pp.par.cov=pp.par.cov,chunk.size=chunk.size)})
	#calculate log pointwise predictive density
	cat("calculating wAIC score...\n\n\n")
		lpd <- calculate.lpd(chunk.lnls,n.iter)
	#calculate effective number of parameters
	pwaic <- calculate.pwaic(chunk.lnls)
	elpd <- lpd - pwaic
	waic <- -2 * elpd
	waic.list <- list("chunk.lnls" = chunk.lnls,"lpd" = lpd,"pwaic" = pwaic,"elpd" = elpd,"waic" = waic)
	return(waic.list)
}

#'@export
compare.x.K <- function(dir,chunk.size){ #FIX TO MAKE FLEXIBLE W/R/T WHICH K DIRS ARE IN DIR
	setwd(dir)
	K.dirs <- list.dirs(full.names=TRUE,recursive=FALSE)
	waics <- vector("list",length=length(K.dirs))
	for(k in 1:length(K.dirs)){
		cat("calculating WAIC for K =",k,"\n\n")
		setwd(K.dirs[k])
		load(list.files(pattern="freq.data"))
		load(list.files(pattern="data.block"))
		load(list.files(pattern="geoStr.results"))
		waics[[k]] <- calculate.waic(freqs = freq.data$freqs,
									 data.block = data.block,
									 geoStr.results = geoStr.results[[1]], #FIX TO MAKE FLEXIBLE W/R/T NO. CHAINS
									 chunk.size=chunk.size,
									 samples=(geoStr.results[[1]]$posterior$n.iter/2 + 1) : geoStr.results[[1]]$posterior$n.iter)
		setwd("..")
	}
	return(waics)
}

#'@export
x.validation <- function(test.pct,n.reps,K,freqs,geoDist,coords,prefix,n.iter,burnin){
	x.val <- lapply(1:n.reps,
					function(i){
						x.validation.rep(rep.no = i,
										 test.pct,
										 K,
										 freqs,
										 geoDist,
										 coords,
										 prefix,
										 n.iter,
										 burnin)
					})
	names(x.val) <- paste0("rep_",1:n.reps)
	return(x.val)
}

#'@export
x.validation.rep <- function(rep.no,test.pct,K,freqs,geoDist,coords,prefix,n.iter,burnin){
	freqs <- drop.invars(freqs)
	train.loci <- sample(1:ncol(freqs),ncol(freqs)*(1-test.pct))
	test.loci <- c(1:ncol(freqs))[!(1:ncol(freqs)) %in% train.loci]
	train.data <- freqs[,train.loci]
	test.data <- freqs[,test.loci]
		save(train.data,file=paste0(prefix,"_rep",rep.no,"_training.dataset.Robj"))
		save(test.data,file=paste0(prefix,"_rep",rep.no,"_testing.dataset.Robj"))
	training.runs.sp <- lapply(K,function(k){
								geoStructure(spatial = TRUE,
											 K = k,
											 freqs = train.data,
											 geoDist = geoDist,
											 coords = coords,
											 prefix = paste0(prefix,"_sp_","rep",rep.no,"K=",k),
											 n.iter = n.iter,
											 burnin = burnin)
						})
	training.runs.nsp <- lapply(K,function(k){
								geoStructure(spatial = FALSE,
											 K = k,
											 freqs = train.data,
											 geoDist = geoDist,
											 coords = coords,
											 prefix = paste0(prefix,"_nsp_","rep",rep.no,"K=",k),
											 n.iter = n.iter,
											 burnin = burnin)
						})
	training.runs <- list("sp" = training.runs.sp,
						  "nsp" = training.runs.nsp)
	save(training.runs,file=paste0(prefix,"rep",rep.no,"_training.runs.Robj"))
	test.lnl.sp <- lapply(training.runs.sp,
						function(x){
							fit.to.test(test.data,x[[1]])
						})
	names(test.lnl.sp) <- K						
	test.lnl.nsp <- lapply(training.runs.nsp,
						function(x){
							fit.to.test(test.data,x[[1]])
						})
	names(test.lnl.nsp) <- K
	test.lnl <- list("sp" = test.lnl.sp,
					 "nsp" = test.lnl.nsp)
	save(test.lnl,file=paste0(prefix,"rep",rep.no,"_test.lnl.Robj"))
	return(test.lnl)
}

#'@export
fit.to.test <- function(test.data,geoStr.results){
	freq.data <- process.freq.data(test.data)
	pp.par.cov <- post.process.par.cov(geoStr.results,
										samples = (1 + geoStr.results$posterior$n.iter/2):
															geoStr.results$posterior$n.iter)
	test.lnl <- lapply(pp.par.cov,
						function(x){
							log.likelihood(freq.data$obsCov,x$inv,x$log.det,freq.data$n.loci)
				})
	return(test.lnl)
}

#'@export
get.seq.max <- function(n.TI.steps){
	seq.max <- 0
	my.max <- 0
	while(my.max < 1){
		my.max <- max(cumsum(1 - 1/(1+seq(0,seq.max,length.out=n.TI.steps))))
		seq.max <- seq.max + 1e-4	
	}
	return(seq.max-1e-4)
}

#'@export
calculate.mean.lnL <- function(geoStr.results,data.block){
	pp.par.cov <- post.process.par.cov(geoStr.results,
										samples = (1 + geoStr.results$posterior$n.iter/2):
													geoStr.results$posterior$n.iter)
	lnLs <- lapply(pp.par.cov,function(x){log.likelihood(data.block$obsCov,x$inv,x$log.det,data.block$L)})
	mean.lnL <- mean(unlist(lnLs))
	return(mean.lnL)
}

#'@export
calculate.model.evidence <- function(temps,mean.lnLs){
	model.evidence <- 0
	for(i in 1:(length(temps)-1)){
		model.evidence <- model.evidence + 0.5 * (mean.lnLs[[i+1]] + mean.lnLs[[i]])*(temps[i+1] - temps[i])
	}
	return(model.evidence)
}

#'@export
nice.rounding <- function(x){
	dig <- 1
	while(any(table(round(x,dig)) > 1)){
		dig <- dig + 1
	}
	return(dig)
}

#'@export
estimate.model.evidence <- function(n.TI.steps,spatial=TRUE,K,freqs,geoDist=NULL,coords,prefix,n.iter=1e3){
	seq.max <- get.seq.max(n.TI.steps)
	temps <- pmin(cumsum(1 - 1 / (1 + seq(0,seq.max,length.out=n.TI.steps))),1)
	
	freq.data <- process.freq.data(freqs)
	data.block <- make.data.block(K,freq.data,coords,spatial,geoDist,temp=NULL)
	geoStr.runs <- lapply(temps,
							function(B){
								geoStructure(spatial,
											 K,
											 freqs,
											 geoDist = geoDist,
											 temp = B,
											 coords = coords,
											 prefix = paste0("temp_",round(B,nice.rounding(temps)),"_",prefix),
											 n.chains = 1,
											 n.iter = n.iter,
											 burnin = 0)
							})
	mean.lnLs <- lapply(geoStr.runs,
						 function(x){
							calculate.mean.lnL(x[[1]],data.block)
						 })
	model.evidence <- calculate.model.evidence(temps,mean.lnLs)
	me.list <- list("temps" = temps,
					"mean.lnLs" = mean.lnLs,
					"data.block" = data.block,
					"model.evidence" = model.evidence)
	save(me.list,file="model.evidence.list.Robj")
	return(model.evidence)
}


