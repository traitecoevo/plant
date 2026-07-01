PACKAGE := $(shell grep '^Package:' DESCRIPTION | sed -E 's/^Package:[[:space:]]+//')
RSCRIPT = Rscript --no-init-file

all: compile

rebuild: clean RcppR6 full_compile roxygen

compile:
	Rscript -e 'pkgbuild::compile_dll(compile_attributes = FALSE, debug=FALSE)' 

debug: RcppR6
	Rscript -e 'pkgbuild::compile_dll(debug=TRUE)' \ 
	make roxygen

# compared to compile, also generates src/RcppExports.cpp, R/RcppExports.R 
full_compile:
	Rscript -e 'pkgbuild::compile_dll(debug=FALSE)' 

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

test: all
	Rscript -e 'library(methods); devtools::test()'

# AD tests self-skip under load_all (pkgload DLL). Run against the INSTALLED
# package: library(plant) loads the installed .so so the skip-gate passes;
# test_dir (not test_local) avoids re-loading via load_all; the namespace-parented
# env exposes unexported internals; parallel off so callr workers inherit the load.
test-ad: install
	TESTTHAT_PARALLEL=false Rscript -e 'library(plant); e <- new.env(parent = asNamespace("plant")); testthat::test_dir("tests/testthat", filter = "ad", env = e)'

benchmark:
	Rscript scripts/benchmark.R

install:
	R CMD INSTALL .

build:
	R CMD build .

check: build
	R CMD check --no-manual `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -f `ls -1tr ${PACKAGE}*gz | tail -n1`
	@rm -rf ${PACKAGE}.Rcheck

clean:
	rm -f src/*.o src/*.so src/*.o.tmp

.PHONY: all compile doc clean test test-ad attributes roxygen install build check
