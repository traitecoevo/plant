TARGETS = $(subst md,pdf,$(shell ls *.md))

all: $(TARGETS)

%.pdf: %.md
	pandoc $<  --include-in-header=include.tex -V linkcolor:black -o output/$@

clean:
	rm -f output/*.pdf
	
