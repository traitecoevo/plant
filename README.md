# TREE model

[![Build Status](http://acropora.bio.mq.edu.au:8080//github.com/richfitz/tree/status.svg?branch=master)](http://acropora.bio.mq.edu.au:8080//github.com/richfitz/tree)

 - TRait Eco-Evolutionary model?
 - TRait Evolutionary Ecology model?
 - Trait Ecosystem Ecology model?
 - ...perhaps some other acronym expansion?
 
# Installation

The entire project directory is an R package (albeit one that is not
valid, because it is not yet fully documented).  The only non-standard
requirement is the Rcpp package and gsl installed where the compilers
can find it.  Then

    R CMD INSTALL tree

will install the package.

There is a Makefile (`GNUmakefile`, to avoid `R CMD SHLIB` during the
compile phase) that is useful for checking compilation without doing a
full install.

# Tests

After installation, tests can be run by running

    make test

from the directory `tree/` or from `inst/tests` will run the tests.
In the latter directory, running

	make test_valgrind

will run tests under valgrind, which can turn up memory errors.

Alternatively, from within R, run

    library(testthat)
	library(tree)
	test_package("tree")

# Documentation

Currently being cobbled together.  Running `make` in the `doc/`
directory will compile some notes to PDF.
