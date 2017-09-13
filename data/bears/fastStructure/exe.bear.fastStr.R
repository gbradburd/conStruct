setwd("/Applications/STRUCTURE/fastStructure")


for(k in 2:7){
	infile <- "/Users/bradburd/Dropbox/conStruct/data/bears/fastStructure/bears"
	outfile <- paste0("/Users/bradburd/Dropbox/conStruct/data/bears/fastStructure/bears_K",k)
	call <- paste0("python structure.py -K ",k," --input=",infile," --output=",outfile," --prior=logistic --format=str")
	system(call)
}

