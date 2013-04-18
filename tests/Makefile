RUN_TESTS="library(testthat); test_dir('.')"

test:
	Rscript -e $(RUN_TESTS)

# Basic command line
VALGRIND=valgrind
# File descriptor for output when doing suppression generation
VALGRIND_FD=4
# Used to toggle full leak checking
VALGRIND_LEAK=--leak-check=full
# Loads the suppression file if available
VALGRIND_OPTS=--suppressions=valgrind.supp  --dsymutil=yes

test_valgrind: valgrind.supp
	R -d '${VALGRIND} ${VALGRIND_OPTS}' \
		--slave --vanilla -e $(RUN_TESTS)

test_valgrind_leak: valgrind.supp
	R -d '${VALGRIND} ${VALGRIND_LEAK} ${VALGRIND_OPTS}' \
		--slave --vanilla -e $(RUN_TESTS)

# Generate the suppressions file 'valgrind.supp'
# 
# This just loads R and the packages that we use; any error there
# cannot be our fault, so we can safely ignore them.
# Options for doing suppression generation
VALGRIND_SUPP_CMD="library(testthat); library(Rcpp); dyn.load('test.so')"
VALGRIND_SUPP=--gen-suppressions=all --log-fd=${VALGRIND_FD}
valgrind.supp:
	echo "void test() {}" > test.c
	R CMD SHLIB test.c
	R -d '${VALGRIND} ${VALGRIND_SUPP}'  \
		--slave --vanilla -e ${VALGRIND_SUPP_CMD}     \
		${VALGRIND_FD}> valgrind.out
	./trim_valgrind.py
	rm -f test.c test.o test.so

clean:
	rm -f test.c test.o test.so

.PHONY: clean test test_valgrind test_valgrind_leak
