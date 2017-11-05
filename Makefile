PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

all: compile

compile: RcppR6 
	Rscript -e 'devtools::compile_dll()' \ 
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
	(cd docs; remake install_vignettes)

staticdocs: vignettes
	@mkdir -p inst/staticdocs
	Rscript -e "library(methods); staticdocs::build_site()"
	rm -f vignettes/*.html
	@rmdir inst/staticdocs
website: staticdocs
	./update_web.sh

.PHONY: all compile doc clean test install vignettes
