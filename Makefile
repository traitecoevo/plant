all:
	cd src; make

clean:
	cd src; make clean

RUN_TESTS="library(testthat); test_dir('.')"
test: all
	cd tests; Rscript -e $(RUN_TESTS)

# Basic command line
VALGRIND=valgrind --dsymutil=yes --suppressions=valgrind-saved.supp
# File descriptor for output when doing suppression generation
VALGRIND_FD=4
# Options for doing suppression generation
VALGRIND_SUPP=--gen-suppressions=all --log-fd=${VALGRIND_FD}
# Uncomment/comment to toggle full leak checking
# VALGRIND_LEAK=--leak-check=full

# Loads the suppression file if available
ifeq ($(wildcard valgrind.supp),) 
    VALGRIND_OPTS=--suppressions=valgrind.supp
endif

test_valgrind: all
	cd tests; R -d '${VALGRIND} ${VALGRIND_LEAK} ${VALGRIND_OPTS}' \
		--slave --vanilla -e $(RUN_TESTS)

# Generate the suppressions file 'valgrind.supp'
# 
# Could possibly improve on this by first doing a load of Rcpp and
# testthat and saving all errors (as they're not my problem).
test_valgrind_suppress: all
	cd tests; R -d '${VALGRIND} ${VALGRIND_LEAK} ${VALGRIND_SUPP}' \
		--slave --vanilla -e $(RUN_TESTS) \
		${VALGRIND_FD}> valgrind.out
	cd tests; ./trim_valgrind.py

.PHONY: clean test test_valgrind test_valgrind_suppress
