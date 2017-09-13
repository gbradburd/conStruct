.PHONY: all, clean

all: conStruct.pdf

conStruct.pdf : reference.library.tex.bib

clean: 
	rm problem_statement-1_0.pdf problem_statement-1.asy problem_statement-1.pre problem_statement-1.tex problem_statement.aux problem_statement.fdb_latexmk problem_statement.fls problem_statement.log problem_statement.pdf problem_statement.pre

%.pdf : %.tex %.bbl
	while ( pdflatex $<;  grep -q "Rerun to get" $*.log ) do true ; done

%.aux : %.tex
	-pdflatex $<

%.bbl : %.aux
	bibtex $<

%.html : %.md
	Rscript -e "templater::render_template(md.file='$<', output='$@')"

%.svg : %.pdf
	inkscape $< --export-plain-svg=$@

%.png : %.pdf
	convert -density 300 $< -flatten $@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

