RUN_TESTS="library(methods); library(testthat); test_dir('.')"

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
# This just loads R and the package that we use; any error there
# cannot be our fault (and there are lots reported on OS X due to the
# false positive problem there), so we can safely ignore them.
VALGRIND_SUPP_CMD="library(testthat); library(tree)"
VALGRIND_SUPP=--gen-suppressions=all --log-fd=${VALGRIND_FD}
valgrind.supp:
	R -d '${VALGRIND} ${VALGRIND_SUPP}'  \
		--slave --vanilla -e ${VALGRIND_SUPP_CMD}     \
		${VALGRIND_FD}> valgrind.out
	./trim_valgrind.py

.PHONY: test test_valgrind test_valgrind_leak
