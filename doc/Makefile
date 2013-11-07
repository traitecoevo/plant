TARGETS = $(subst md,pdf,$(shell ls *.md)) Growth_model_pars.pdf

all: $(TARGETS)


ML-challenges.pdf: figs/fitnesslandscape/plot.pdf

Growth_model_pars.md: Growth_model_pars.csv
	Rscript -e "library(pander, quietly =TRUE); x<- read.csv('$<', stringsAsFactors=FALSE, header = TRUE); cat(pandoc.table.return(x, styl='multiline', split.tables= 200, split.cells = 50, justify = 'left'), file='$@')"

figs/fitnesslandscape/plot.pdf: figs/fitnesslandscape/fitnesslandscape.r
	Rscript $<

%.pdf: %.md
	pandoc $<  --include-in-header=include.tex -V linkcolor:black -o $@

clean:
	rm -f *.pdf
