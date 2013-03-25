all:
	cd src; make

clean:
	cd src; make clean

RUN_TESTS="library(testthat); test_dir('.')"
test: all
	cd tests; Rscript -e $(RUN_TESTS)

# Basic command line
VALGRIND=valgrind --dsymutil=yes
# File descriptor for output when doing suppression generation
VALGRIND_FD=4
# Options for doing suppression generation
VALGRIND_SUPP=--gen-suppressions=all --log-fd=${VALGRIND_FD}
# Used to toggle full leak checking
VALGRIND_LEAK=--leak-check=full
# Loads the suppression file if available
VALGRIND_OPTS=--suppressions=valgrind-saved.supp

test_valgrind: all tests/valgrind-saved.supp
	cd tests; R -d '${VALGRIND} ${VALGRIND_OPTS}' \
		--slave --vanilla -e $(RUN_TESTS)

test_valgrind_leak: all tests/valgrind-saved.supp
	cd tests; R -d '${VALGRIND} ${VALGRIND_LEAK} ${VALGRIND_OPTS}' \
		--slave --vanilla -e $(RUN_TESTS)

# Generate the suppressions file 'valgrind.supp'
# 
# This just loads R and the packages that we use; any error there
# cannot be our fault, so we can safely ignore them.
tests/valgrind-saved.supp: all
	cd tests; R -d '${VALGRIND} ${VALGRIND_LEAK} ${VALGRIND_SUPP}'  \
		--slave --vanilla -e "library(testthat); library(Rcpp)" \
		${VALGRIND_FD}> valgrind.out
	cd tests; ./trim_valgrind.py
	cat tests/valgrind.supp >> tests/valgrind-saved.supp

.PHONY: clean test test_valgrind test_valgrind_leak
