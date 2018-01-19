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
	#recover()
	n.snps <- (ncol(ped)-6)/2
	map <- matrix(NA,nrow=n.snps,ncol=3)
	#CHR
	map[,1] <- rep(1,n.snps)
	#RS
	map[,2] <- paste0("rs",1:n.snps)
	#BP
	map[,3] <- seq(1,n.snps*1e3,length.out=n.snps)
	return(map)
}

ped_poplar <- convert.str.to.ped(str.file="~/Dropbox/conStruct/data/poplars/fastStructure/poplars.str")
map_poplar <- make.map.file(ped_poplar)
	write.table(ped_poplar,file="~/Dropbox/conStruct/data/poplars/admixture/poplar.ped",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(map_poplar,file="~/Dropbox/conStruct/data/poplars/admixture/poplar.map",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
call <- "plink --noweb --file poplar --make-bed --out poplar --map3 --allow-no-sex"
system(call)
file.remove("~/Dropbox/conStruct/data/poplars/admixture/poplar.ped")
file.remove("~/Dropbox/conStruct/data/poplars/admixture/poplar.map")
file.remove("~/Dropbox/conStruct/data/poplars/admixture/poplar.log")
file.remove("~/Dropbox/conStruct/data/poplars/admixture/poplar.nosex")