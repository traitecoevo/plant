PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

all: compile

rebuild: clean RcppR6 full_compile roxygen

compile:
	Rscript -e 'pkgbuild::compile_dll(compile_attributes = FALSE)' 

# compared to compile, also generates src/RcppExports.cpp, R/RcppExports.R 
full_compile:
	Rscript -e 'pkgbuild::compile_dll()' 

# generates 
RcppR6:
	Rscript -e "library(methods); RcppR6::RcppR6()"

# generates src/RcppExports.cpp, R/RcppExports.R from anything with Rcpp::export. 
attributes:
	Rscript -e "Rcpp::compileAttributes()"

# generates documentation
roxygen:
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"

test:
	Rscript -e 'library(methods); devtools::test()'

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

slow_vignettes:
	(cd inst/slow_vignettes; ln -sfn ../../vignettes vignettes; remake update_vignettes)

vignettes:
	Rscript -e "devtools::build_vignettes()"

website: vignettes
	Rscript -e "pkgdown::build_site()" /
	open "inst/website/index.html"

.PHONY: all compile doc clean test attributes roxygen install build check vignettes website push_website
