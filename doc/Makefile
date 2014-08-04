TARGETS = $(subst md,pdf,$(shell ls *.md))

all: $(TARGETS)

%.pdf: %.md include.tex
	pandoc $< --toc --template=include.tex --variable fontsize=12pt --latex-engine=xelatex -o $@

clean:
	rm -f *.pdf
