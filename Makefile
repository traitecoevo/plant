PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')

all:
# This is probably The Right Way to do this:
#	Rscript -e 'devtools::compile_dll()'
# but will overzealosuly compile everything.
	cd src; R CMD SHLIB *.cpp -o ${PACKAGE}.so

test:
	Rscript -e 'library(methods); devtools::test()'

RcppR6:
	Rscript -e "library(methods); RcppR6::RcppR6()"

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
