setwd("/Applications/STRUCTURE/fastStructure")

K <- 2
infile <- "/Users/bradburd/Dropbox/conStruct/sims/structure/datasets/simK3/simK3"
outfile <- paste0("/Users/bradburd/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK",K)
call <- paste0("python structure.py -K ",K," --input=",infile," --output=",outfile," --prior=logistic --format=str")
system(call)

K <- 3
infile <- "/Users/bradburd/Dropbox/conStruct/sims/structure/datasets/simK3/simK3"
outfile <- paste0("/Users/bradburd/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK",K)
call <- paste0("python structure.py -K ",K," --input=",infile," --output=",outfile," --prior=logistic --format=str")
system(call)

K <- 4
infile <- "/Users/bradburd/Dropbox/conStruct/sims/structure/datasets/simK3/simK3"
outfile <- paste0("/Users/bradburd/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK",K)
call <- paste0("python structure.py -K ",K," --input=",infile," --output=",outfile," --prior=logistic --format=str")
system(call)