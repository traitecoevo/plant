all: EBT.pdf details.pdf

%.pdf: %.md
	pandoc $<  --include-in-header=include.tex -V linkcolor:black -o $@

clean:
	rm -f EBT.pdf
