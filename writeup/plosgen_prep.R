################################################################
################################################################
#	ridiculous script that prepares 
#	paper for plos genetics submission
################################################################
################################################################

# Call this script from terminal

setwd("~/Dropbox/conStruct/writeup")

# make sure flags are:
#	submissionversiontrue
#	includesupplementtrue
#	floatsatendfalse
#	includefigstrue
call <- c("pdflatex conStruct.tex")
system(call)
system(call)

# split out suppmat figs
call <- c("pdfseparate -f 31 -l 58 conStruct.pdf S%d.pdf")
system(call)

# rename suppmat figs
supp.files <- list.files(pattern="^S")
for(i in 31:58){
	file.rename(paste0("S",i,".pdf"),paste0("S",i-30,"_fig.pdf"))
}

# read in text, switch flags, and recompile
text <- scan("conStruct.tex",what=character(),sep="\n")
new.text <- text
new.text[grepl("\\includesupplementtrue",text)] <- "%\\includesupplementtrue"
new.text[grepl("\\includesupplementfalse",text)] <- "\\includesupplementfalse"
new.text[grepl("\\\\floatsatendfalse",text)] <- "%\\floatsatendfalse"
new.text[grepl("\\\\floatsatendtrue",text)] <- "\\floatsatendtrue"
	cat(new.text,file="conStruct.tex",sep="\n")

call <- c("pdflatex conStruct.tex")
system(call)
system(call)
system(call)

# split out main text figs and tables
call <- c("pdfseparate -f 26 -l 35 conStruct.pdf Fig%d.pdf")
system(call)

call <- c("pdfseparate -f 37 -l 37 conStruct.pdf Table%d.pdf")
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
new.text <- text
new.text[grepl("\\includesupplementtrue",text)] <- "\\includesupplementtrue"
new.text[grepl("\\includesupplementfalse",text)] <- "%\\includesupplementfalse"
new.text[grepl("\\\\floatsatendfalse",text)] <- "\\floatsatendfalse"
new.text[grepl("\\\\floatsatendtrue",text)] <- "%\\floatsatendtrue"
new.text[grepl("\\\\includefigsfalse",text)] <- "\\includefigsfalse"
new.text[grepl("\\\\includefigstrue",text)] <- "%\\includefigstrue"
	cat(new.text,file="conStruct.tex",sep="\n")

call <- c("pdflatex conStruct.tex")
system(call)
system(call)

# split out the pre-suppmat pages
call <- c("pdfseparate -f 1 -l 25 conStruct.pdf Page%d.pdf")
system(call)

# unite the main text into a single doc for submission
call <- paste0("pdfunite ",c(paste0("Page",1:25,".pdf",collapse=" "))," conStruct.pdf")
system(call)

file.remove(paste0("Page",1:25,".pdf"))