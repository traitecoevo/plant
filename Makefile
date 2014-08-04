all:
	make -C src

doc:
	make -C doc

clean:
	make -C doc clean
	make -C src clean

test:
	make -C inst/tests test

test-regression:
	make -C tests/regression

## This one takes a long time.
scripts:
	make -C scripts

# I dislike devtools' use of the Collate field (which causes problems
# for rapidly adding new code) so I'm disabling it this way:
DEVTOOLS_DOCUMENT=devtools::document(roclets=c('namespace', 'rd'))
roxygen: all
	@mkdir -p man
	Rscript -e "library(methods); ${DEVTOOLS_DOCUMENT}"

RSRC = $(shell find ./R)

install:
	R CMD INSTALL .

.PHONY: all doc clean test install
