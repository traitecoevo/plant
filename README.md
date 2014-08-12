# The TRait Ecology and Evolution (TREE) model

[![Build Status](http://acropora.bio.mq.edu.au:8080//github.com/richfitz/tree/status.svg?branch=master)](http://acropora.bio.mq.edu.au:8080//github.com/richfitz/tree)

# Installation

The entire project directory is an R package (albeit one that is not
valid, because it is not yet fully documented).  The only non-standard
requirement is that the gnu scientific library (gsl) is installed where the compilers
can find it.  In addition you'll need some other R packages (see list [here](https://github.com/richfitz/tree/issues/83)). Then run

    make install

to compile and install the package.

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
directory will compile some notes to PDF. These require installation
of the [pandoc package](http://johnmacfarlane.net/pandoc/installing.html).
