PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

all: compile

compile: RcppR6
	Rscript -e 'pkgbuild::compile_dll()' \ 
	make roxygen

test:
	Rscript -e 'library(methods); devtools::test()'

RcppR6:
	Rscript -e "library(methods); RcppR6::RcppR6()"

attributes:
	Rscript -e "Rcpp::compileAttributes()"

roxygen:
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"

benchmark:
	Rscript scripts/benchmark.R

install:
	R CMD INSTALL .

build:
	R CMD build --no-build-vignettes .

check: build
	R CMD check --no-build-vignettes --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck

clean:
	rm -f src/*.o src/*.so

vignettes:
	(cd inst/docs; ln -sfn ../../vignettes vignettes; remake install_vignettes)

website: vignettes
	Rscript -e "pkgdown::build_site()" /
	open "inst/website/index.html"

.PHONY: all compile doc clean test attributes roxygen install build check vignettes website push_website
