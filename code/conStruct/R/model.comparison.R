#'@export
x.validation <- function(test.pct,n.reps,K,freqs,geoDist,coords,prefix,n.iter){
	x.val <- lapply(1:n.reps,
					function(i){
						x.validation.rep(rep.no = i,
										 test.pct,
										 K,
										 freqs,
										 geoDist,
										 coords,
										 prefix,
										 n.iter)
					})
	names(x.val) <- paste0("rep_",1:n.reps)
	return(x.val)
}

#'@export
x.validation.rep <- function(rep.no,test.pct,K,freqs,geoDist,coords,prefix,n.iter,make.figs=FALSE,save.files=FALSE){
	freqs <- drop.invars(freqs)
	train.loci <- sample(1:ncol(freqs),ncol(freqs)*(1-test.pct))
	test.loci <- c(1:ncol(freqs))[!(1:ncol(freqs)) %in% train.loci]
	train.data <- freqs[,train.loci]
	test.data <- freqs[,test.loci]
		save(train.data,file=paste0(prefix,"_rep",rep.no,"_training.dataset.Robj"))
		save(test.data,file=paste0(prefix,"_rep",rep.no,"_testing.dataset.Robj"))
	training.runs.sp <- lapply(K,function(k){
								conStruct(spatial = TRUE,
											 K = k,
											 freqs = train.data,
											 geoDist = geoDist,
											 coords = coords,
											 prefix = paste0(prefix,"_sp_","rep",rep.no,"K",k),
											 n.iter = n.iter,
											 make.figs = make.figs,
											 save.files = save.files)
						})
	training.runs.nsp <- lapply(K,function(k){
								conStruct(spatial = FALSE,
											 K = k,
											 freqs = train.data,
											 geoDist = geoDist,
											 coords = coords,
											 prefix = paste0(prefix,"_nsp_","rep",rep.no,"K",k),
											 n.iter = n.iter,
											 make.figs = make.figs,
											 save.files = save.files)
						})
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
match.clusters.x.runs <- function(csr1,csr2,csr1.order=NULL){
	cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	K1 <- ncol(csr1$MAP$admix.proportions)
		if(!is.null(csr1.order)){
			csr1$MAP$admix.proportions <- csr1$MAP$admix.proportions[,csr1.order]
		}
		K1.cols <- cluster.colors[1:K1]
	K2 <- ncol(csr2$MAP$admix.proportions)
		K2.cols <- numeric(K2)
	k.combn <- expand.grid(1:K1,1:K2)
	clst.sims <- unlist(lapply(1:nrow(k.combn),
						function(n){
							measure.frob.similarity(csr1$MAP$admix.proportions[,k.combn[n,1],drop=FALSE],
													csr2$MAP$admix.proportions[,k.combn[n,2],drop=FALSE],
													K=1)
							}))
	while(length(which(K2.cols == 0)) > (K2-K1)){
		tmp.max <- which.max(rank(clst.sims,na.last=FALSE))
		csr2.match <- k.combn[tmp.max,2]
		csr1.match <- k.combn[tmp.max,1]
		K2.cols[csr2.match] <- K1.cols[csr1.match]
		clst.sims[which(k.combn[,1]==csr1.match)] <- NA
		clst.sims[which(k.combn[,2]==csr2.match)] <- NA
	}
	if(K2 > K1){
		K2.cols[which(K2.cols==0)] <- cluster.colors[(K1+1):K2]
	}
	K2.clst.order <- unique(match(cluster.colors,K2.cols)[which(!is.na(match(cluster.colors,K2.cols)))])
	clst2.info <- list("cols" = K2.cols,
					   "clst.order" = K2.clst.order)
	return(clst2.info)
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


determine.log.shift <- function(chunk.lnls){
	if(diff(range(unlist(chunk.lnls))) > 700){
		message("the difference between the min and max lnLs may be inducing underflow")
	}
	return(max(unlist(chunk.lnls)))
}


shift.chunk.lnls <- function(chunk.lnls,A,n.iter){
	shift.chunk.lnls <- log(sum(exp(chunk.lnls-A))) + A - log(n.iter)
	return(shift.chunk.lnls)
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