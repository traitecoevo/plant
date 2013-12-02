TARGETS = $(subst md,pdf,$(shell ls *.md)) Growth_model_pars.pdf

all: $(TARGETS)


ML-challenges.pdf: figs/fitnesslandscape/plot.pdf

Growth_model_pars.md: Growth_model_pars.csv
	Rscript -e "library(pander, quietly =TRUE); x<- read.csv('$<', stringsAsFactors=FALSE, header = TRUE); cat(pandoc.table.return(x, styl='multiline', split.tables= 200, split.cells = 50, justify = 'left'), file='$@')"

figs/fitnesslandscape/plot.pdf: figs/fitnesslandscape/fitnesslandscape.r
	Rscript $<

%.pdf: %.md include.tex
	pandoc $< --toc --template=include.tex --variable fontsize=12pt --latex-engine=xelatex -o $@

clean:
	rm -f *.pdf
