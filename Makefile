all:
	make -C doc
	make -C src

clean:
	make -C doc clean
	make -C src clean

test: all
	make -C inst/tests test

document: all
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"
