TARGETS = output/$(subst md,pdf,$(shell ls *.md))

all: $(TARGETS)

ML-challenges.pdf: figs/fitnesslandscape/plot.pdf

figs/fitnesslandscape/plot.pdf: figs/fitnesslandscape/fitnesslandscape.r
	Rscript $<

%.pdf: %.md
	pandoc $<  --include-in-header=include.tex -V linkcolor:black -o output/$@

clean:
	rm -f output/*.pdf
	
