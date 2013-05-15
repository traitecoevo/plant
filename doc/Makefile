all: EBT.pdf

EBT.pdf: EBT.md
	pandoc $<  --include-in-header=include.tex -V linkcolor:black -o $@

clean:
	rm -f EBT.pdf
