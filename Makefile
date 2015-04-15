all:
	make -C src

doc:
	make -C doc

clean:
	make -C doc clean
	make -C src clean

test:
	make -C inst/tests test

test-regression:
	make -C tests/regression

## This one takes a long time.
scripts:
	make -C scripts

attributes:
	Rscript -e "Rcpp::compileAttributes()"

roxygen:
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"

install:
	R CMD INSTALL .

.PHONY: all doc clean test install
