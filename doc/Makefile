TARGETS = $(subst md,pdf,$(shell ls *.md)) $(subst Rmd,pdf,$(shell ls *.Rmd))

all: $(TARGETS)

%.md: %.Rmd
	Rscript -e "knitr::knit('$<')"

%.pdf: %.md include.tex the-american-naturalist.csl
	pandoc $< --csl=the-american-naturalist.csl --toc --filter pandoc-citeproc --bibliography=refs.bib --template=include.tex --variable fontsize=12pt --latex-engine=xelatex -o $@

the-american-naturalist.csl:
	curl -k https://raw.githubusercontent.com/citation-style-language/styles/master/the-american-naturalist.csl > $@

clean:
	rm -f *.pdf
