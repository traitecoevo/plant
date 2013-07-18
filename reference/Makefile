TARGETS = single
# This is where the repository falster-traitdiversity is located on my
# (RGF) computer.  Run with
#   PATH_EVOLVE = /path/to/falster-traitdiversity make
# to override.  No need to edit this file.
PATH_EVOLVE ?= ~/Documents/Projects/veg/falster-traitdiversity

all: ${TARGETS}

single:
	PATH_EVOLVE=${PATH_EVOLVE} Rscript -e\
		"library(testthat); test_file('reference.R')"

clean:
	rm -rf ref-single
