################################################################
################################################################
#	ridiculous script that prepares 
#	paper for plos genetics submission
################################################################
################################################################

# Call this script from terminal

setwd("~/Dropbox/conStruct/writeup")

# first, make sure flags are:
#	submissionversiontrue
#	includesupplementtrue
#	floatsatendfalse
#	includefigstrue
text <- readLines("conStruct.tex") 
text[which(text == "\\submissionversiontrue")] <- "\\submissionversiontrue"
text[which(text == "%\\submissionversiontrue")] <- "\\submissionversiontrue"
text[which(text == "\\submissionversionfalse")] <- "%\\submissionversionfalse"
text[which(text == "%\\submissionversionfalse")] <- "%\\submissionversionfalse"

text[which(text == "\\includesupplementtrue")] <- "\\includesupplementtrue"
text[which(text == "%\\includesupplementtrue")] <- "\\includesupplementtrue"
text[which(text == "\\includesupplementfalse")] <- "%\\includesupplementfalse"
text[which(text == "%\\includesupplementfalse")] <- "%\\includesupplementfalse"

text[which(text == "\\floatsatendtrue")] <- "%\\floatsatendtrue"
text[which(text == "%\\floatsatendtrue")] <- "%\\floatsatendtrue"
text[which(text == "\\floatsatendfalse")] <- "\\floatsatendfalse"
text[which(text == "%\\floatsatendfalse")] <- "\\floatsatendfalse"

text[which(text == "\\includefigstrue")] <- "\\includefigstrue"
text[which(text == "%\\includefigstrue")] <- "\\includefigstrue"
text[which(text == "\\includefigsfalse")] <- "%\\includefigsfalse"
text[which(text == "%\\includefigsfalse")] <- "%\\includefigsfalse"
writeLines(text,"conStruct.tex")

# compile
call <- c("pdflatex conStruct.tex")
system(call)
system(call)

# run git-latexdiff to make a diff file 
#	to the initial submission 55053bbd120f68ed66e199bb39ecd30d3686dc19
call <- c("git-latexdiff --main conStruct.tex 55053bbd120f68ed66e199bb39ecd30d3686dc19 --bibtex --ignore-makefile -o sub1_vs_sub2_diff.pdf")
system(call)

# split out suppmat figs
call <- c("pdfseparate -f 32 -l 62 conStruct.pdf S%d.pdf")
system(call)

# rename suppmat figs
supp.files <- list.files(pattern="^S")
for(i in 32:62){
	file.rename(paste0("S",i,".pdf"),paste0("S",i-31,"_fig.pdf"))
}

# read in text, switch flags, and recompile
text <- readLines("conStruct.tex")
text[which(text == "\\includesupplementtrue")] <- "%\\includesupplementtrue"
text[which(text == "%\\includesupplementfalse")] <- "\\includesupplementfalse"
text[which(text == "\\floatsatendfalse")] <- "%\\floatsatendfalse"
text[which(text == "%\\floatsatendtrue")] <- "\\floatsatendtrue"
writeLines(text,"conStruct.tex")

call <- c("pdflatex conStruct.tex")
system(call)
system(call)
system(call)

# split out main text figs and tables
call <- c("pdfseparate -f 28 -l 37 conStruct.pdf Fig%d.pdf")
system(call)

call <- c("pdfseparate -f 39 -l 39 conStruct.pdf Table%d.pdf")
system(call)

# rename main text figs and tables
fig.files <- list.files(pattern="^Fig")
for(i in 1:length(fig.files)){
	file.rename(fig.files[i],paste0("Fig",i,".pdf"))
}

table.files <- list.files(pattern="^Table")
for(i in 1:length(table.files)){
	file.rename(table.files[i],paste0("Table",i,".pdf"))
}

# prepare file to make submission pdf
#	no floats at end, no figs (but still captions),
#	and suppmat, so the internal references don't break
text <- readLines("conStruct.tex")
text[which(text == "%\\includesupplementtrue")] <- "\\includesupplementtrue"
text[which(text == "\\includesupplementfalse")] <- "%\\includesupplementfalse"
text[which(text == "%\\floatsatendfalse")] <- "\\floatsatendfalse"
text[which(text == "\\floatsatendtrue")] <- "%\\floatsatendtrue"
text[which(text == "%\\includefigsfalse")] <- "\\includefigsfalse"
text[which(text == "\\includefigstrue")] <- "%\\includefigstrue"
writeLines(text,"conStruct.tex")

call <- c("pdflatex conStruct.tex")
system(call)
system(call)

# split out the pre-suppmat pages
call <- c("pdfseparate -f 1 -l 26 conStruct.pdf Page%d.pdf")
system(call)

# unite the main text into a single doc for submission
call <- paste0("pdfunite ",c(paste0("Page",1:26,".pdf",collapse=" "))," conStruct.pdf")
system(call)

file.remove(paste0("Page",1:26,".pdf"))