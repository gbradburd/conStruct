.PHONY: all, clean, epsfigs

all: conStruct.pdf

conStruct.pdf : reference.library.tex.bib

clean: 
	-rm *.aux *.log *.lof *.lot *.fff *.ttt *.out *.bbl *.blg

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

%.eps : %.pdf
	inkscape --without-gui --export-eps=$@ $<

epsfigs : 
	-mkdir -p eps_figs
	for x in  $$(grep includegr conStruct.tex | sed -e 's/.*{//' | sed -e 's/}//g'); do y=$$(echo $$x | sed -e 's/pdf/eps/'); make $$y; mv $$y eps_figs; done
