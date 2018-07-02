convert.str.to.ped <- function(str.file){
	str <- read.table(str.file,header=FALSE,stringsAsFactors=FALSE,row.names=1)
	n.loci <- ncol(str)-1
	n.haps <- nrow(str)
	hap.IDs <- row.names(str)
	pop.IDs <- as.numeric(unlist(lapply(lapply(strsplit(hap.IDs,"_"),"[[",2),function(x){strsplit(x,"\\.")[[1]][[1]]}))[seq(1,n.haps,by=2)])
	mand.cols <- matrix(0,nrow=n.haps/2,ncol=6)
	#family
	mand.cols[,1] <- pop.IDs
	#ind ID
	mand.cols[,2] <- unlist(lapply(unique(pop.IDs),function(x){1:length(which(pop.IDs==x))}))
	genos <- matrix(NA,nrow=n.haps/2,ncol=n.loci*2)
	genos[,seq(1,n.loci*2,by=2)] <- as.matrix(str[seq(1,n.haps,by=2),2:ncol(str)])
	genos[,seq(2,n.loci*2,by=2)] <- as.matrix(str[seq(2,n.haps,by=2),2:ncol(str)])
	genos <- genos + 1
	ped <- cbind(mand.cols, genos)
	return(ped)
}

make.map.file <- function(ped){
	n.snps <- (ncol(ped)-6)/2
	map <- matrix(NA,nrow=n.snps,ncol=4)
	#CHR
	map[,1] <- rep(1,n.snps)
	#RS
	map[,2] <- paste0("rs",1:n.snps)
	#morgans
	map[,3] <- rep(0,n.snps)
	#BP
	map[,4] <- seq(1,n.snps*1e3,length.out=n.snps)
	return(map)
}

ped_K1 <- convert.str.to.ped(str.file="simK1.str")
map_K1 <- make.map.file(ped_K1)
	write.table(ped_K1,file="simK1/simK1.ped",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(map_K1,file="simK1/simK1.map",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
setwd("simK1")
call <- "plink --noweb --file simK1 --out simK1"
system(call)
file.remove("simK1.ped")
file.remove("simK1.map")
file.remove("simK1.log")
file.remove("simK1.nosex")

setwd("../simK2")
ped_K2 <- convert.str.to.ped(str.file="simK2.str")
map_K2 <- make.map.file(ped_K2)
	write.table(ped_K2,file="simK2/simK2.ped",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(map_K2,file="simK2/simK2.map",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
call <- "plink --file simK2 --out simK2"
system(call)
file.remove("simK2.ped")
file.remove("simK2.map")
file.remove("simK2.log")
file.remove("simK2.nosex")


setwd("../simK3")
ped_K3 <- convert.str.to.ped(str.file="simK3.str")
map_K3 <- make.map.file(ped_K3)
	write.table(ped_K3,file="simK3/simK3.ped",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(map_K3,file="simK3/simK3.map",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
call <- "plink --file simK3 --out simK3"
system(call)
file.remove("simK3.ped")
file.remove("simK3.map")
file.remove("simK3.log")
file.remove("simK3.nosex")