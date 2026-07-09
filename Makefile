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

benchmark:
	Rscript scripts/benchmark.R

# Run the TF24 hydraulic scenario gateway and write a scorecard RDS.
scenarios:
	Rscript scripts/run_scenario_gateway.R

# Re-bless the gateway baseline: regenerate the scorecard the opt-in gateway
# test (tests/testthat/test-scenario-gateway.R) diffs against. Run only when a
# changed outcome is intended.
bless-scenarios:
	SCENARIO_OUT=tests/testthat/test_data/scenario_baseline.rds \
		Rscript scripts/run_scenario_gateway.R

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

.PHONY: all compile doc clean test attributes roxygen install build check benchmark scenarios bless-scenarios
