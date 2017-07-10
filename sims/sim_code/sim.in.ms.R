##################################################################
##################################################################
##	run ms on a lattice
##################################################################
##################################################################

require(conStruct)

within.cluster.m <- function(sample.inds,migration.index,pairwise.migration.matrix){
	#recover()
	migration.rate.vector <- c()
	for(i in 1:nrow(migration.index)){
		migration.rate.vector <- c(migration.rate.vector,
									sprintf("-m %s %s %s",
										sample.inds[migration.index[i,1]],
										sample.inds[migration.index[i,2]],
										pairwise.migration.matrix[migration.index[i,1],migration.index[i,2]]))
	}
	return(migration.rate.vector)
}

within.cluster.merge <- function(sample.inds,shallow.split){
	merge.vector <- unlist(lapply(sample.inds[2]:sample.inds[length(sample.inds)],
						function(i){
							sprintf("-ej %s %s %s",
									shallow.split,
									i,
									sample.inds[1])}))
	return(merge.vector)
}

rewrite.admix.props <- function(w){
	K <- length(w)
	new.w <- w
	for(k in 2:(K-1)){
		adjustment <- (1-sum(w[1:(k-1)]))
		if(adjustment == 0){
			new.w[k] <- w[k]
		} else {
			new.w[k] <- w[k]/(1-sum(w[1:(k-1)]))
		}
	}
	new.w[new.w > 1] <- 1
	return(new.w[1:(K-1)])
}

write.admixture.event.calls <- function(admix.list,n.pops){
	#recover()
	N <- nrow(admix.list$w)
	K <- ncol(admix.list$w)
	admixture.call <- c()
	source.ticker <- 1
	n.time.ticker <- 1e-5/K
	k.time.ticker <- 1e-5/(10*K)
	for(n in 1:N){
		if(K > 2){
			new.w <- rewrite.admix.props(admix.list$w[n,])
		} else if (K == 2){
			new.w <- admix.list$w[n,1]
		}
		for(i in 1:length(new.w)){
			admixture.call <- c(admixture.call,
								sprintf("-es %s %s %s",
											admix.list$time + n*n.time.ticker + i*k.time.ticker,
											admix.list$targets[[n]][i],
											new.w[i]),
								sprintf("-ej %s %s %s",
											admix.list$time + n*n.time.ticker + i*k.time.ticker + 1e-10,
											n.pops + source.ticker,
											admix.list$sources[[n]][i]))
			source.ticker <- source.ticker + 1
		}
	}
	return(admixture.call)
}

btwn.cluster.merge <- function(K,n.pops,deep.split){
	if(K > 1){
		merge.vector <- unlist(lapply(2:K,
						function(k){
							sprintf("-ej %s %s %s",
								deep.split,
								1+(k-1)*n.pops,
								1)
						}))
	} else {
		merge.vector <- c()
	}
	return(merge.vector)
}

write.migration.rates <- function(K,coords,migration.rate,shallow.split,deep.split,admix.list=NULL){
	# recover()
	pop.dist <- fields::rdist(coords)
	pop.dist[which(pop.dist > sqrt(2) | pop.dist < 1)] <- Inf
	pairwise.migration.matrix <- lapply(1:K,function(k){migration.rate[k]/pop.dist})
	migration.index <- lapply(pairwise.migration.matrix,function(x){which(x != 0,arr.ind=TRUE)})
	n.pops <- nrow(coords)
	sample.inds <- lapply(1:K,function(k){1:n.pops + (k-1) * n.pops})
	migration.rate.vector <- c()
	for(k in 1:K){
		migration.rate.vector <- c(migration.rate.vector,
									within.cluster.m(sample.inds[[k]],migration.index[[k]],pairwise.migration.matrix[[k]]))
	}
	for(k in 1:K){
		migration.rate.vector <- c(migration.rate.vector,
									within.cluster.merge(sample.inds[[k]],shallow.split[k]))
	}
	migration.rate.vector <- c(migration.rate.vector,
									btwn.cluster.merge(K,n.pops,deep.split))
	if(!is.null(admix.list)){
		migration.rate.vector <- c(migration.rate.vector,
									write.admixture.event.calls(admix.list,n.pops*K))
	}
	return(migration.rate.vector)
}

# code cannibalized from Dan Denison that reads ms output into R
read.ms.haplotype.matrices <- function(nsam, ndraws, ms.output.file) {
    txt <- scan(file=ms.output.file, what=character(0), sep="\n", quiet=TRUE)
    h <- list()
    ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
    marker <- grep("segsites", txt)
    stopifnot(length(marker) == ndraws)
    ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
    segsites <- sapply(strsplit(txt[marker], split=":"), function(vec) as.integer(vec[2]))
    for(draw in seq(along=marker)) {
        if(!(draw %% 100)) cat(draw, " ")
        if(segsites[draw] > 0) {
            haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
            haplotypes <- strsplit(haplotypes, split="")
            h[[draw]] <- sapply(haplotypes, function(el) c(as.integer(el)))
            ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
            if(segsites[draw] == 1) h[[draw]] <- as.matrix(h[[draw]])
            ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
            else h[[draw]] <- t(h[[draw]])
        }
        else h[[draw]] <- matrix(nrow=nsam, ncol=0)
        stopifnot(all(dim(h[[draw]]) == c(nsam, segsites[draw])))  
    }
    cat("\n")
    h
}

write.lattice <- function(coords,K,n.chromo,subsample=NULL){
	#recover()
	if(is.null(subsample)){
		to.sample <- rep(1,nrow(coords))
	} else {
		to.sample <- numeric(nrow(coords))
		to.sample[subsample] <- 1
	}
	to.sample <- c(to.sample,rep(0,(K-1)*nrow(coords))) * n.chromo
	to.sample <- paste0(to.sample,collapse=" ")
	return(to.sample)
}

# function that actually calls ms with the call put together by all these R functions
ms <- function(coords,K,n.chromo,ms.params,admix.list=NULL,subsample=NULL){
		ms.output.file <- "ms_output"
		random.seeds <- c(sample(1:100000,3,replace=TRUE))
		n.pops <- nrow(coords)
		sampled.pops <- ifelse(is.null(subsample),n.pops,length(subsample))
		call <- paste(
					sprintf(
						"/Applications/ms.folder/msdir/ms %s 1 -t %s -s 1 -I %s %s 0.0-m %s -seeds %s %s %s",
							sampled.pops*n.chromo,
							ms.params$theta,
							K*n.pops,
							write.lattice(coords,K,n.chromo,subsample),
							paste0(write.migration.rates(K,
														 coords,
														 ms.params$migration.rate,
														 ms.params$shallow.split,
														 ms.params$deep.split,
														 admix.list),
														 collapse=" "),
							random.seeds[1],
							random.seeds[2],
							random.seeds[3]
					),
					"-T", ">", ms.output.file
				)
		cat(call,file="ms_call.txt")
		system(call)
		read.ms.haplotype.matrices(nsam=sampled.pops*n.chromo,ndraws=1,ms.output.file=ms.output.file)
}

# generate times and migration rates in ms units
generate.ms.command.line.values <- function(diploid.population.size,locus.size,per.bp.mu,migration.fraction,deep.split,shallow.split){
	ms.command.line.values <- vector("list",length=4)
		names(ms.command.line.values) <- c("theta","deep.split","shallow.split","migration.rate")
			ms.command.line.values$theta <- 4*diploid.population.size*per.bp.mu*locus.size
			ms.command.line.values$migration.rate <- 4*diploid.population.size*migration.fraction
			ms.command.line.values$deep.split <- deep.split/(4*diploid.population.size)
			ms.command.line.values$shallow.split <- shallow.split/(4*diploid.population.size)
			# ms.command.line.values$time <- 4*diploid.population.size*generations.ago
			#admixture on time scale more recent than 1/(4Nm k) won't spread
	return(ms.command.line.values)
}

check.param <- function(param,K){
	if(length(param != K)){
		param <- rep(param,K)
	}
	return(param)
}

divvy.space <- function(radius,centroid,K){
	points <- matrix(NA,K,2)
	angles <- 3*pi/4 + seq(0,2*pi,by=2*pi/K)[1:K]
	for(k in 1:K){
		points[k,] <- centroid + radius * c(cos(angles[k]),sin(angles[k]))
	}
	return(points)
}

adjudicate.borders <- function(fdr,K,in.K.dist){
	#recover()
	w <- numeric(K)
	if(any(fdr < in.K.dist)){
		w[which.min(fdr)] <- 1
	} else {
		if(runif(1) > 0.5){
			w[which.min(fdr)] <- 1
		} else {
			dcps <- (1/fdr)^2
			if(K > 2){
				dcps[which.max(fdr)] <- 0
			}
			dcps <- dcps/(2*max(dcps))
			w <- gtools::rdirichlet(1,dcps)
		}
	}
	return(w)
}

get.nearest.neighbs <- function(coords){
	nns <- lapply(1:nrow(coords),
					function(n){
						n.dists <- fields::rdist(coords[n,,drop=FALSE],coords) ;
						n.dists[which(n.dists==0)] <- NA ;
						which(n.dists == min(n.dists,na.rm=TRUE))
					})
	return(nns)
}

tidy.w.sim <- function(w,N,coords,focal.dist.ratios){
	#recover()
	nns <- get.nearest.neighbs(coords)
	for(n in 1:N){
		if(!any(w[n,] == 1)){
			focal.max <- which.max(w[n,])
			neighb.max <- apply(w[nns[[n]],,drop=FALSE],1,which.max)
			ticker <- 0
			while(!any(neighb.max==focal.max) & ticker < 50){
				fdr <- unlist(lapply(focal.dist.ratios,"[[",n))
				dcps <- (1/fdr)^2
				dcps[which.max(fdr)] <- 0
				dcps <- dcps/(2*max(dcps))
				w[n,] <- gtools::rdirichlet(1,dcps)
				focal.max <- which.max(w[n,])
				ticker <- ticker + 1
			}
		}
	}
	return(w)
}

sim.admix.props <- function(N,K,coords,subsample=NULL){
	#recover()
	w <- matrix(1,N,K)
	if(is.null(subsample)){
		pops.to.check <- 1:nrow(coords)
	} else {
		pops.to.check <- subsample
	}
	while(!any(apply(w[pops.to.check,],1,function(x){any((x^2)!=x)}))){
		if(is.null(subsample)){
			radius <- sqrt(2*(diff(range(coords[,1]))/2)^2)
		} else {
			radius <- sqrt(2*(diff(range(coords[subsample,1]))/2)^2)
		}
		K.foci <- divvy.space(radius = radius,
						      centroid=colMeans(coords),
						  	  K=K)	
		dist.to.foci <- lapply(1:nrow(K.foci),
							function(k){
								fields::rdist(K.foci[k,,drop=FALSE],coords)
							})
		focal.dist.ratios <- lapply(1:K,function(k){c(dist.to.foci[[k]]/Reduce("+",dist.to.foci))})
		in.K.dist <- quantile(unlist(focal.dist.ratios),0.1)
		for(n in 1:N){
			w[n,] <- adjudicate.borders(unlist(lapply(focal.dist.ratios,"[[",n)),K,in.K.dist)
		}
		w <- tidy.w.sim(w,N,coords,focal.dist.ratios)
	}
	return(w)
}

make.admix.list <- function(N,K,coords,w,time){
	admix.list <- list("w" = w,
					   "N" = N,
					   "K" = K,
					   "time" = time,
					   "sources" = lapply(1:N,function(n){n + N*(seq(1,K-1,by=1))}),
					   "targets" = lapply(1:N,function(n){n + N*(seq(0,K-2,by=1))}))
	return(admix.list)
}

subsample.admix.list <- function(admix.list,subsample){
	admix.list <- list("w" = admix.list$w[subsample,],
					   "N" = length(subsample),
					   "K" = admix.list$K,
					   "time" = admix.list$time)
	return(admix.list)
}

# make dataset for use by, e.g., conStruct
generate.conStruct.dataset <- function(n.loci,K,coords,n.chromo,ms.params,admix.list=NULL,subsample=NULL){
#	recover()
	#Allele Counts & Sample sizes
		ms.params$migration.rate <- check.param(ms.params$migration.rate,K)
		ms.params$shallow.split <- check.param(ms.params$shallow.split,K)
		data.matrix <- do.call(cbind,
							replicate(n.loci,
								ms(coords,K,n.chromo,ms.params,admix.list,subsample)))
		n.pops <- ifelse(is.null(subsample),nrow(coords),length(subsample))
		population.membership <- unlist(lapply(1:n.pops,function(n){rep(n,n.chromo)}))
		allele.counts <- matrix(0,nrow=n.pops,ncol=n.loci)
		for(i in 1:n.pops){
			allele.counts[i,] <- colSums(data.matrix[which(population.membership==i),,drop=FALSE])
		}
		sample.sizes <- matrix(n.chromo,nrow=n.pops,ncol=n.loci)
	#Return sim output
		if(!is.null(subsample)){
			admix.list <- subsample.admix.list(admix.list,subsample)
			coords <- coords[subsample,]
		}
		N <- nrow(coords)		
		sim.dataset <- list("ms.params" = ms.params,
							"admix.list" = admix.list,
							"freq.data" = list("counts" = allele.counts,
											   "sample.sizes" = sample.sizes,
											   "freqs" = allele.counts/sample.sizes),
							"coords" = coords,
							"N" = N,
							"K" = K)
	return(sim.dataset)
}

admix.pie.plot <- function(w, coords, cluster.colors, radii = 2.7, add = FALSE, x.lim = NULL, y.lim = NULL) {
	require(caroline)
	K <- ncol(w)
	N <- nrow(w)
	cluster.names <- paste0("cluster_", 1:K)
    sample.names <- paste0("sample_", 1:N)
    color.tab <- nv(c(cluster.colors[1:K]), cluster.names)
    pie.list <- lapply(1:N, function(i){nv(w[i, ], cluster.names)})
    names(pie.list) <- sample.names
    if (add) {
        par(new = TRUE)
    }
    else {
        par(mar = c(2, 2, 2, 2))
    }
    if (is.null(x.lim)) {
        x.lim <- c(min(coords[, 1]) - 1, max(coords[, 1]) + 1)
    }
    if (is.null(y.lim)) {
        y.lim <- c(min(coords[, 2]) - 1, max(coords[, 2]) + 1)
    }
    pies(pie.list, x0 = coords[, 1], y0 = coords[, 
        2], color.table = color.tab, border = "black", radii = radii, 
        xlab = "", ylab = "", main = "", lty = 1, density = NULL, 
        xlim = x.lim, ylim = y.lim)
    box(lwd = 2)
    return(invisible(0))
}

