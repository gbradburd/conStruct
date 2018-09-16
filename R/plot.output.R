#' Make output plots
#'
#' \code{make.all.the.plots} makes figures from the output from a 
#' 	conStruct analysis.
#'
#' This function takes the output from a conStruct analysis and 
#' generates a number of plots for visualizing results and 
#' diagnosing MCMC performance.
#'
#' @param conStruct.results The list output by a 
#'			\code{conStruct} run.
#' @param data.block A \code{data.block} list saved during a 
#'			\code{conStruct} run.
#' @param prefix A character vector to be prepended to all figures.
#' @param layer.colors A \code{vector} of colors to be used in 
#'			plotting results for different layers. Users must 
#'			specify one color per layer.  If \code{NULL}, plots 
#'			will use a pre-specified vector of colors.
#' @return This function has only invisible return values.
#'
#'	@details This function produces a variety of plots that can be 
#'	useful for visualizing results or diagnosing MCMC performance. 
#'  The plots made are by no means an exhaustive, and users are 
#' 	encouraged to make further plots, or customize these plots as they 
#'	see fit.  For each plot, one file is generated for each MCMC chain 
#'	(specified with the \code{n.chains} argument in the function 
#'	\code{conStruct}. The plots generated (as .pdf files) are:
#'	\itemize{
#'		\item Structure plot - STRUCTURE-style plot, where each sample 
#'			is represented as a stacked bar plot, and the length of the 
#'			bar plot segments of each color represent that sample's 
#'			admixture proportion in that layer. Described further 
#'			in the help page for \code{make.structure.plot}.
#'		\item Admixture pie plot - A map of samples in which each sample's 
#'				location is denoted with a pie chart, and the proportion 
#'				of a pie chart of each color represents that sample's 
#'				admixture in each layer. Described further in the help 
#'				page for \code{make.admix.pie.plot}
#'		\item model.fit.CIs - A plot of the sample allelic covariance 
#'			shown with the 95\% credible interval of the parametric 
#'			covariance for each entry in the matrix.
#'		\item layer.covariances - A plot of the layer-specific 
#'				covariances overlain unto the sample allelic covariance.
#'		\item Trace plots - Plots of parameter values over the MCMC.
#'		\itemize{
#'			\item lpd - A plot of the log posterior probability over the MCMC.
#'			\item nuggets - A plot of estimates of the nugget parameters 
#'				over the MCMC.
#'			\item gamma - A plot of estimates of the gamma parameter 
#'				over the MCMC.
#'			\item layer.cov.params - Plots of estimates of the 
#'				layer-specific parameters over the MCMC.
#'			\item admix.props - A plot of estimates of the admixture proportions 
#'				over the MCMC.
#'		}
#'	}
#' 
#'@export
make.all.the.plots <- function(conStruct.results,data.block,prefix,layer.colors=NULL){
	if(!any(grepl("chain",names(conStruct.results)))){
		stop("\nyou must specify conStruct results across all chains\ni.e. from conStruct.results rather than conStruct.results[[1]]\n\n")
	}
	if(!is.null(layer.colors)){
		if(length(layer.colors!=data.block$K)){
			stop("\nyou must specify one color per layer\n\n")
		}
	} else {
		layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
		if(data.block$K > 10){
			stop("\nyou has specified more layers than there are default colors.\n you must specify your own layer.colors")
		}
	}
	lapply(1:length(conStruct.results),function(i){
		make.all.chain.plots(conStruct.results[[i]],chain.no=i,data.block,prefix,layer.colors)
	})
	return(invisible("made chain plots!"))
}

#' Make STRUCTURE output plot
#'
#' \code{make.structure.plot} makes a STRUCTURE-style plot from the output from a 
#' 	conStruct analysis.
#'
#' This function takes the output from a conStruct analysis and 
#' makes a STRUCTURE-style plot, where each sample 
#' is represented as a stacked bar plot, and the length of the 
#' bar plot segments of each color represent that sample's 
#' admixture proportion in that layer.
#'
#' @param admix.proportions A \code{matrix} of admixture proportions, 
#'			with one row per sample and one column per layer.
#' @param mar A \code{vector} of plotting margins passed to \code{par}.
#'			Default is \code{c(2,4,2,2)}, which tends to look good.
#' @param sample.order A \code{vector} giving the order in which sample 
#'			admixture proportions are to be plotted, left to right.  If 
#'			\code{NULL}, samples are plotted in the order they occur in 
#'			\code{admix.proportions}.
#' @param layer.order A \code{vector} giving the order in which layers 
#'			are plotted, bottom to top. If \code{NULL}, layers are plotted 
#'			in the order they occur in \code{admix.proportions}.
#' @param sample.names Vector of names to be plotted under each sample's 
#'			admixture proportion bar plot. The index of a sample's name 
#'			should be the same as the index of the sample's row in 
#'			\code{admix.proportions}. If \code{NULL}, no names 
#'			are printed.
#' @param sort.by An \code{integer} giving the column index of the \code{admix.proportions} 
#'			matrix to be used in determining sample plotting order.  If specified, 
#'			samples are plotted from left to right in increasing order of their 
#'			membership in that layer.  If \code{NULL}, samples are plotted 
#'			in the order they occur in \code{admix.proportions}.
#' @param layer.colors A \code{vector} of colors to be used in plotting 
#'			results for different layers. Users must specify one 
#'			color per layer.  If \code{NULL}, the plot will use 
#'			a pre-specified vector of colors.
#' @return This function has only invisible return values.
#'
#' @examples
#' \dontshow{
#'		admix.props <- matrix(c(0.086,0.000,0.500,0.505,0.099,0.052,0.024,0.007,0.800,0.000,0.216,0.744,0.917,0.199,0.469,0.000,0.783,0.298,0.329,0.446,0.000,0.000,0.637,0.903,0.000,0.000,0.000,0.012,0.021,0.000,0.000,0.089,0.000,0.554,0.002,0.000,0.000,0.095,0.020,0.001,0.001,0.011,0.000,0.200,0.000,0.060,0.053,0.082,0.036,0.013,0.000,0.062,0.169,0.137,0.029,0.001,0.000,0.178,0.079,0.000,0.999,1.000,0.988,0.979,0.975,1.000,0.744,0.984,0.435,0.998,0.914,1.000,0.405,0.475,0.900,0.947,0.965,0.993,0.000,1.000,0.725,0.203,0.000,0.765,0.518,1.000,0.154,0.533,0.534,0.525,0.999,1.000,0.185,0.018,1.000,0.001,0.000,0.000,0.000,0.025,0.000,0.167,0.016,0.012,0.000),ncol=3)
#'	}	
#' # make STRUCTURE-style plot
#' 	make.structure.plot(admix.proportions = admix.props)
#'
#' # make STRUCTURE-style plot, sorted by membership in layer 1
#'	make.structure.plot(admix.proportions = admix.props,sort.by=1) 
#' 
#'@export
make.structure.plot <- function(admix.proportions,mar=c(2,4,2,2),sample.order=NULL,layer.order=NULL,sample.names=NULL,sort.by=NULL,layer.colors=NULL){
	if(class(admix.proportions)!="matrix"){
		stop("\nyou must specify a matrix of admixture proportions\n")
	}
	K <- ncol(admix.proportions)
	N <- nrow(admix.proportions)
	graphics::par(mar=mar)
	if(is.null(layer.order)){
		layer.order <- seq(1:K)
	}
	if(is.null(sample.order)){
		sample.order <- seq(1:N)
	}
	if(!is.null(sort.by)){
		if(sort.by > K){
			stop("\nyou must specify a layer that exists in your data\n")
		}
		sample.order <- order(admix.proportions[,sort.by])
	}
	if(is.null(layer.colors)){
		layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	} else {
		if(K > length(layer.colors)){
			stop("\nyou must specify one color per layer\n")
		}
	}
	use.colors <- layer.colors[1:K][layer.order]
	graphics::plot(0,xlim=c(0,N),ylim=c(0,1),type='n',ylab="admixture",xlab="",xaxt='n')
	plotting.admix.props <- apply(cbind(0,admix.proportions[,layer.order]),1,cumsum)
	lapply(1:K,function(i){
		make.structure.polygon.layer(plotting.admix.props,i,use.colors,sample.order)
	})
	if(!is.null(sample.names)){
		graphics::axis(side=1,at=seq(1:N)-0.5,labels=sample.names[sample.order],cex.axis=0.5,las=2)
	}
	return(invisible("plotted"))
}

#' Make admixture pie plot
#'
#' \code{make.structure.plot} makes a map of pie plots showing admixture 
#' proportions across layers.
#'
#' This function takes the output from a conStruct analysis and 
#' makes a map of pie plots showing admixture proportions across layers, 
#' where each sample is represented as a pie chart, and the proportion of 
#' the pie of each color represent that sample's 
#' admixture proportion in that layer.
#'
#' @param admix.proportions A \code{matrix} of admixture proportions, 
#'			with one row per sample and one column per layer.
#' @param coords \code{matrix} of sample coordinates, with one row 
#'			per sample and two columns giving (respectively) the X 
#'			and Y plotting coordinates.
#' @param layer.colors A \code{vector} of colors to be used in 
#'			plotting results for different layers. Users must 
#'			specify one color per layer.  If \code{NULL}, the plot 
#'			will use a pre-specified vector of colors.
#' @param radii A \code{vector} of numeric values giving the radii to be 
#'			used in plotting admixture pie plots. If the number of values 
#'			specified is smaller than the number of samples, radii values 
#'			will be recycled across samples. The default is 2.7.
#' @param add A \code{logical} value indicating whether to add the pie plots 
#'			to an existing plot.  Default is \code{FALSE}.
#' @param x.lim A \code{vector} giving the x limits of the plot. The default
#'			value is \code{NULL}, which indicates that the range of values 
#'			given in the first column of \code{coords} should be used.
#' @param y.lim A \code{vector} giving the y limits of the plot. The default
#'			value is \code{NULL}, which indicates that the range of values 
#'			given in the second column of \code{coords} should be used.
#' @param mar A \code{vector} giving the number of lines of margin specified 
#'		for the four sides of the plotting window (passed to \code{par}). 
#'		Default value, which is only used if \code{add=FALSE}, is 
#'		\code{c(2,2,2,2)}.
#' @return This function has only invisible return values.
#' @examples
#' \dontshow{
#'		admix.props <- matrix(c(0.086,0.000,0.500,0.505,0.099,0.052,0.024,0.007,0.800,0.000,0.216,0.744,0.917,0.199,0.469,0.000,0.783,0.298,0.329,0.446,0.000,0.000,0.637,0.903,0.000,0.000,0.000,0.012,0.021,0.000,0.000,0.089,0.000,0.554,0.002,0.000,0.000,0.095,0.020,0.001,0.001,0.011,0.000,0.200,0.000,0.060,0.053,0.082,0.036,0.013,0.000,0.062,0.169,0.137,0.029,0.001,0.000,0.178,0.079,0.000,0.999,1.000,0.988,0.979,0.975,1.000,0.744,0.984,0.435,0.998,0.914,1.000,0.405,0.475,0.900,0.947,0.965,0.993,0.000,1.000,0.725,0.203,0.000,0.765,0.518,1.000,0.154,0.533,0.534,0.525,0.999,1.000,0.185,0.018,1.000,0.001,0.000,0.000,0.000,0.025,0.000,0.167,0.016,0.012,0.000),ncol=3)
#'		coords <- matrix(c(-126.38,-125.23,-126.97,-128.54,-126.95,-121.71,-126.79,-123.38,-137.88,-125.82,-122.94,-130.73,-123.08,-122.84,-128.58,-124.82,-129.75,-122.25,-122.32,-129.10,-125.28,-123.98,-133.35,-131.74,-124.16,-146.35,-94.63,-149.02,-111.50,-126.67,-133.77,-118.63,-115.78,-113.42,-135.33,52.40,49.84,54.66,54.65,51.69,49.44,52.82,50.05,59.52,51.34,45.81,56.81,44.71,50.24,54.14,51.04,56.68,52.98,54.04,55.34,50.64,50.23,58.76,57.30,50.54,64.90,56.35,63.87,56.92,65.23,68.38,54.75,60.80,50.82,60.70),ncol=2)
#'	}	
#' # make admixture pie plot
#'	make.admix.pie.plot(admix.proportions = admix.props,coords = coords)
#' 
#'@export
make.admix.pie.plot <- function(admix.proportions,coords,layer.colors=NULL,radii=2.7,add=FALSE,x.lim=NULL,y.lim=NULL,mar=c(2,2,2,2)){
	if(class(admix.proportions)!="matrix"){
		stop("\nyou must specify a matrix of admixture proportions\n")
	}
	if(is.null(coords)){
		stop("\nyou must specify sampling coordinates\n")
	} else {
		if(nrow(coords) != nrow(admix.proportions)){
			stop("\nyou must specify one set of coordinates for each row of the admixture proportion matrix\n")
		}
	}
	K <- ncol(admix.proportions)
	N <- nrow(admix.proportions)	
	if(is.null(layer.colors)){
		layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	} else {
		if(K > length(layer.colors)){
			stop("\nyou must specify one color per layer\n")
		}
	}
	layer.names <- paste0("layer_",1:K)
	sample.names <- paste0("sample_",1:N)
	color.tab <- caroline::nv(c(layer.colors[1:K]),layer.names)
	pie.list <- lapply(1:N,function(i){caroline::nv(admix.proportions[i,],layer.names)})
	names(pie.list) <- sample.names
	if(add){
		graphics::par(new=TRUE)
	} else {
		graphics::par(mar=mar)
	}
	if(is.null(x.lim)){
		x.lim <- c(min(coords[,1]) - 1, max(coords[,1]) + 1)
	}
	if(is.null(y.lim)){
		y.lim <- c(min(coords[,2]) - 1, max(coords[,2]) + 1)
	}
	suppressWarnings(
		caroline::pies(pie.list,x0=coords[,1],y0=coords[,2],
					   color.table=color.tab,border="black",radii=radii,
					   xlab="",ylab="",main="",lty=1,density=NULL,
					   xlim = x.lim, ylim = y.lim)
	)
	return(invisible(0))
}

#' Compare two conStruct runs
#'
#' \code{compare.two.runs} makes figures comparing the output 
#' 	from two conStruct analyses.
#'
#' This function takes the outputs from two conStruct analyses and 
#' generates a number of plots for comparing results and 
#' diagnosing MCMC performance.
#'
#' @param conStruct.results1 The list output by a 
#'			\code{conStruct} run.
#' @param data.block1 A \code{data.block} list saved during a 
#'			\code{conStruct} run.
#' @param conStruct.results2 The list output by a second
#'			\code{conStruct} run.
#' @param data.block2 A \code{data.block} list saved during a 
#'			second \code{conStruct} run.
#' @param prefix A character vector to be prepended to all figures.
#' @param layer.colors A \code{vector} of colors to be used in 
#'			plotting results for different layers. Users must 
#'			specify one color per layer.  If \code{NULL}, plots 
#'			will use a pre-specified vector of colors.
#' @return This function has only invisible return values.
#'
#'	@details This function produces a variety of plots that can be 
#'	useful for comparing results from two \code{conStruct} analyses.
#'  The runs must have the same number of independent MCMC chains, 
#' 	but may have different values of \code{K}. The spatial and 
#'	nonspatial models can be compared. If the runs were executed 
#'	with different values of \code{K}, the run with the smaller 
#'	value of \code{K} should be specified in the first set of 
#'	arguments (\code{conStruct.results1} and \code{data.block1}).
#'
#'  The plots made are by no means an exhaustive, and users are 
#' 	encouraged to make further plots, or customize these plots as they 
#'	see fit.  For each plot, one file is generated for each MCMC chain 
#'	in each analysis (specified with the \code{n.chains} argument in 
#'	the function \code{conStruct}. For clarity, the layers in the second 
#' 	are matched to those in the first using the function 
#'	\code{match.clusters.x.runs} The plots generated (as .pdf files) are:
#'	\itemize{
#'		\item Structure plot - STRUCTURE-style plot, where each sample 
#'			is represented as a stacked bar plot, and the length of the 
#'			bar plot segments of each color represent that sample's 
#'			admixture proportion in that layer. Described further 
#'			in the help page for \code{make.structure.plot}.
#'		\item Admixture pie plot - A map of samples in which each sample's 
#'				location is denoted with a pie chart, and the proportion 
#'				of a pie chart of each color represents that sample's 
#'				admixture in each layer. Described further in the help 
#'				page for \code{make.admix.pie.plot}
#'		\item model.fit.CIs - A plot of the sample allelic covariance 
#'			shown with the 95\% credible interval of the parametric 
#'			covariance for each entry in the matrix.
#'		\item layer.covariances - A plot of the layer-specific 
#'				covariances overlain unto the sample allelic covariance.
#'		\item Trace plots - Plots of parameter values over the MCMC.
#'		\itemize{
#'			\item lpd - A plot of the log posterior probability over the MCMC.
#'			\item nuggets - A plot of estimates of the nugget parameters 
#'				over the MCMC.
#'			\item gamma - A plot of estimates of the gamma parameter 
#'				over the MCMC.
#'			\item layer.cov.params - Plots of estimates of the 
#'				layer-specific parameters over the MCMC.
#'			\item admix.props - A plot of estimates of the admixture proportions 
#'				over the MCMC.
#'		}
#'	}
#'@export
compare.two.runs <- function(conStruct.results1,data.block1,conStruct.results2,data.block2,prefix,layer.colors=NULL){
	if(length(conStruct.results1) != length(conStruct.results1)){
		stop("\nthe two \"conStruct.results\" objects must be from analyses with the same number of chains\n\n")
	}
	if(data.block1$K > data.block2$K){
		stop("\nyou must specify the run with the smaller value of K in the first set of arguments (\"conStruct.results1\" and \"data.block1\")\n\n")
	}
	if(!any(grepl("chain",names(conStruct.results1)))){
		stop("\nyou must specify conStruct results across all chains\ni.e. from \"conStruct.results\" rather than \"conStruct.results[[1]]\"\n\n")
	}
	if(!is.null(layer.colors)){
		if(length(layer.colors!=data.block1$K)){
			stop("\nyou must specify one color per layer\n\n")
		}
	} else {
		layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
		if(data.block1$K > 10){
			stop("\nyou has specified more layers than there are default colors.\n you must specify your own \"layer.colors\"\n\n")
		}
	}
	lapply(1:length(conStruct.results1),function(i){
		make.all.chain.coplots(conStruct.results1[[i]],conStruct.results2[[i]],chain.no=i,data.block1,data.block2,prefix,layer.colors)
	})
}

plot.lpd <- function(conStruct.results){
	graphics::plot(conStruct.results$posterior$lpd,
			ylab="posterior probability",
			main="Posterior probability",type='l',
			xlab="MCMC iterations")
	return(invisible(0))
}


plot.nuggets <- function(conStruct.results){
	graphics::matplot(conStruct.results$post$nuggets,type='l',
				main="sample nuggets",
				ylab="nugget value",
				xlab="MCMC iterations")
	return(invisible("nuggets"))
}

plot.gamma <- function(conStruct.results){
	graphics::plot(conStruct.results$posterior$gamma,
			ylab="gamma",
			xlab="MCMC iterations",
			main="Gamma",type='l')
	return(invisible(0))
}


get.ylim <- function(layer.params,n.layers,param){
	y.lim <- range(unlist(
				lapply(
					lapply(1:n.layers,
						function(i){
							layer.params[[i]][[param]]
						}),
					function(x){
						range(x)
					})))
	y.lim <- y.lim + c(-0.15*diff(y.lim),0.15*diff(y.lim))
	return(y.lim)
}


plot.layer.param <- function(layer.param,layer.col){
	graphics::points(layer.param,type='l',col=layer.col)
	return(invisible(0))
}


plot.layer.cov.params <- function(data.block,conStruct.results,layer.colors){
	n.layers <- data.block$K
	params <- names(conStruct.results$posterior$layer.params$layer_1)[!names(conStruct.results$posterior$layer.params$layer_1)=="layer.cov"]
	param.ranges <- lapply(params,function(x){get.ylim(conStruct.results$posterior$layer.params,n.layers,x)})
	if(length(params) > 0){
		for(i in 1:length(params)){
			graphics::plot(0,type='n',main=params[i],
				xlab="MCMC iterations",ylab="parameter value",
				ylim=param.ranges[[i]],xlim=c(1,conStruct.results$posterior$n.iter))
			lapply(1:n.layers,function(j){plot.layer.param(conStruct.results$posterior$layer.params[[j]][[params[i]]],layer.colors[j])})
			graphics::legend(x="topright",col= layer.colors[1:n.layers],lty=1,legend=paste0("Layer_",1:n.layers))
		}
	}
	return(invisible(0))
}


plot.admix.props <- function(data.block,conStruct.results,layer.colors,layer.order=NULL){
	n.layers <- data.block$K
	if(is.null(layer.order)){
		layer.order <- 1:n.layers
	}
	graphics::par(mfrow=c(n.layers,1),mar=c(3,3,2,2))
		for(i in 1:n.layers){
			graphics::matplot(conStruct.results$posterior$admix.proportions[,,layer.order[i]],type='l',ylim=c(0,1),
					main=paste0("layer ",i),ylab="admixture proportion",col=layer.colors[i])
		}
	return(invisible(0))
}


get.par.cov.CI <- function(data.block,conStruct.results){
	combns <- gtools::combinations(n=data.block$N,r=2,v=1:data.block$N,repeats.allowed=TRUE)
	CIs <- lapply(1:nrow(combns),
				function(i){
					stats::quantile(conStruct.results$posterior$par.cov[,combns[i,1],combns[i,2]],c(0.025,0.975))
				})
	return(CIs)
}

plot.model.fit.CIs <- function(data.block,conStruct.results){
	cov.range <- range(c(data.block$obsCov,
						conStruct.results$posterior$par.cov))
	graphics::plot(data.block$geoDist,data.block$obsCov,
    	xlab = "geographic distance", 
        ylab = "covariance",
        main="Cov/geoDist",
        ylim = cov.range, type = "n")
	combns <- gtools::combinations(n=data.block$N,r=2,v=1:data.block$N,repeats.allowed=TRUE)
	CIs <- get.par.cov.CI(data.block,conStruct.results)
	lapply(1:nrow(combns),
			function(i){
				graphics::segments(x0 = data.block$geoDist[combns[i,1],combns[i,2]],
						 y0 = CIs[[i]][1],
						 x1 = data.block$geoDist[combns[i,1],combns[i,2]],
						 y1 = CIs[[i]][2],
						 col = grDevices::adjustcolor(1,0.1),
						 lwd=1.5)
			})
	graphics::points(data.block$geoDist,data.block$obsCov,col=2,pch=20,cex=0.8)
	graphics::legend(x="topright",legend=c("observed","95% CI"),pch=c(19,NA),lty=c(NA,1),col=c(2,"gray"))
	return(invisible("plotted"))
}

plot.layer.covariances <- function(data.block,conStruct.results,layer.colors,layer.order=NULL){
	n.layers <- data.block$K
	if(is.null(layer.order)){
		layer.order <- 1:n.layers
	}
	ind.mat <- upper.tri(data.block$geoDist,diag=TRUE)
	order.mat <- order(data.block$geoDist)
	    y.range <- range(c(
	    			unlist(lapply(1:data.block$K,
	    							function(k){
								        conStruct.results$MAP$layer.params[[k]]$layer.cov
					})) + conStruct.results$MAP$gamma, 
					data.block$obsCov))
	graphics::plot(data.block$geoDist[ind.mat],
		 data.block$obsCov[ind.mat],
			xlim=range(data.block$geoDist),ylim=y.range,
			xlab = "geographic distance",
			ylab = "covariance",
			pch=19,col=grDevices::adjustcolor(1,0.7))
		lapply(1:data.block$K, function(k) {
             graphics::lines(data.block$geoDist[order.mat][ind.mat],
             		conStruct.results$MAP$gamma + 
             		conStruct.results$MAP$layer.params[[layer.order[k]]]$layer.cov[order.mat][ind.mat],
                  col = 1,lwd=4.5,lty=1) ; 
             graphics::lines(data.block$geoDist[order.mat][ind.mat],
             		conStruct.results$MAP$gamma + 
             		conStruct.results$MAP$layer.params[[layer.order[k]]]$layer.cov[order.mat][ind.mat],
                  col = layer.colors[k],lwd=4,lty=1)
        })
		graphics::legend(x="topright",col= layer.colors[1:data.block$K],lty=1,
				legend=paste0("Layer_",1:data.block$K),cex=0.7)
	return(invisible("layer covs"))	
}


structure.polygon <- function(plotting.admix.props,i,j,use.colors){
	graphics::polygon(x = c(j-1,j,j,j-1),
				y = c(plotting.admix.props[i,j],
						plotting.admix.props[i,j],
						plotting.admix.props[i+1,j],
						plotting.admix.props[i+1,j]),
				col=use.colors[i])
	return(invisible(j))
}


make.structure.polygon.layer <- function(plotting.admix.props,i,use.colors,sample.order){
		lapply(1:ncol(plotting.admix.props),function(j){
			structure.polygon(plotting.admix.props[,sample.order],i,j,use.colors)
		})
	return(invisible(i))
}

make.all.chain.plots <- function(conStruct.results,chain.no,data.block,prefix,layer.colors){
	grDevices::pdf(file=paste0(prefix,"_trace.plots.chain_",chain.no,".pdf"))
		plot.lpd(conStruct.results)
		plot.nuggets(conStruct.results)
		plot.gamma(conStruct.results)
		plot.layer.cov.params(data.block,conStruct.results,layer.colors)
		if(data.block$K > 1){
			plot.admix.props(data.block,conStruct.results,layer.colors)
		}
	grDevices::dev.off()
	if(!is.null(data.block$geoDist)){
		grDevices::pdf(file=paste0(prefix,"_model.fit.CIs.chain_",chain.no,".pdf"))
			plot.model.fit.CIs(data.block,conStruct.results)
		grDevices::dev.off()
	}
	if(data.block$spatial | (data.block$K > 1 & !is.null(data.block$geoDist))){
		grDevices::pdf(file=paste0(prefix,"_layer.cov.curves.chain_",chain.no,".pdf"),width=5,height=5)
			plot.layer.covariances(data.block,conStruct.results,layer.colors)
		grDevices::dev.off()
	}
	if(data.block$K > 1){
		grDevices::pdf(file=paste0(prefix,"_pie.map.chain_",chain.no,".pdf"),width=6,height=6)	
			make.admix.pie.plot(conStruct.results$MAP$admix.proportions,data.block$coords,layer.colors,radii =2.7,add=FALSE,x.lim=NULL,y.lim=NULL)
		grDevices::dev.off()
		grDevices::pdf(file=paste0(prefix,"_structure.plot.chain_",chain.no,".pdf"),width=10,height=5)
			make.structure.plot(conStruct.results$MAP$admix.proportions,mar=c(2,4,2,2),sample.order=NULL,layer.order=NULL,sample.names=NULL,sort.by=NULL,layer.colors)
		grDevices::dev.off()
	}
	return(invisible("made chain plots!"))
}

make.all.chain.coplots <- function(conStruct.results1,conStruct.results2,chain.no,data.block1,data.block2,prefix,layer.colors){
	match.order <- match.layers.x.runs(conStruct.results1$MAP$admix.proportions,conStruct.results2$MAP$admix.proportions)
	layer.colors1 <- layer.colors
	layer.colors2 <- layer.colors1[match.order]
	grDevices::pdf(file=paste0(prefix,"_trace.plots.chain_",chain.no,".pdf"),width=12,height=6)
		graphics::par(mfrow=c(1,2))
			plot.lpd(conStruct.results1)
			plot.lpd(conStruct.results2)
		graphics::par(mfrow=c(1,2))
			plot.nuggets(conStruct.results1)
			plot.nuggets(conStruct.results2)
		graphics::par(mfrow=c(1,2))
			plot.gamma(conStruct.results1)
			plot.gamma(conStruct.results2)
		graphics::par(mfrow=c(1,2))
			plot.layer.cov.params(data.block1,conStruct.results1,layer.colors1)
			plot.layer.cov.params(data.block2,conStruct.results2,layer.colors2)
		if(data.block1$K > 1 & data.block2$K > 1){
			graphics::par(mfrow=c(1,2))
				plot.admix.props(data.block1,conStruct.results1,layer.colors1)
				plot.admix.props(data.block2,conStruct.results2,layer.colors1,layer.order=match.order)
		}
	grDevices::dev.off()
	if(!is.null(data.block1$geoDist) & !is.null(data.block2$geoDist)){
		grDevices::pdf(file=paste0(prefix,"_model.fit.CIs.chain_",chain.no,".pdf"),width=12,height=6)
			graphics::par(mfrow=c(1,2))
				plot.model.fit.CIs(data.block1,conStruct.results1)
				plot.model.fit.CIs(data.block2,conStruct.results2)
		grDevices::dev.off()
	}
	if(data.block1$spatial & data.block2$spatial | (data.block1$K > 1 & !is.null(data.block1$geoDist) & data.block2$K > 1 & !is.null(data.block2$geoDist))){
		grDevices::pdf(file=paste0(prefix,"_layer.cov.curves.chain_",chain.no,".pdf"),width=10,height=5)
			graphics::par(mfrow=c(1,2))
				plot.layer.covariances(data.block1,conStruct.results1,layer.colors1)
				plot.layer.covariances(data.block2,conStruct.results2,layer.colors1,layer.order=match.order)
		grDevices::dev.off()
	}
	if(data.block1$K > 1 & data.block2$K > 1){
		grDevices::pdf(file=paste0(prefix,"_pie.map.chain_",chain.no,".pdf"),width=12,height=6)	
		graphics::par(mfrow=c(1,2))
			make.admix.pie.plot(conStruct.results1$MAP$admix.proportions,data.block1$coords,layer.colors1,radii =2.7,add=FALSE,x.lim=NULL,y.lim=NULL)
			make.admix.pie.plot(conStruct.results2$MAP$admix.proportions,data.block2$coords,layer.colors1[order(match.order)],radii =2.7,add=FALSE,x.lim=NULL,y.lim=NULL)
		grDevices::dev.off()
		grDevices::pdf(file=paste0(prefix,"_structure.plot.chain_",chain.no,".pdf"),width=10,height=10)
			graphics::par(mfrow=c(2,1))
				make.structure.plot(conStruct.results1$MAP$admix.proportions,mar=c(2,4,2,2),sample.order=NULL,layer.order=NULL,sample.names=NULL,sort.by=NULL,layer.colors1)
				make.structure.plot(conStruct.results2$MAP$admix.proportions,mar=c(2,4,2,2),sample.order=NULL,layer.order=match.order,sample.names=NULL,sort.by=NULL,layer.colors1[order(match.order)])
		grDevices::dev.off()
	}
	return(invisible("made chain plots!"))
}