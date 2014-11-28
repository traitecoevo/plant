PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')

all:
	make -C src

doc:
	make -C doc

clean:
	make -C doc clean
	make -C src clean

test:
	Rscript -e 'library(methods); devtools::test()'

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

check: build
	R CMD check --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck

.PHONY: all doc clean test install
