PKG = $(shell Rscript -e "writeLines(system.file(package='tree'))")

all:
	make -C src

doc:
	make -C doc

clean:
	make -C doc clean
	make -C src clean

test: ${PKG}
	make -C inst/tests test

# I dislike devtools' use of the Collate field (which causes problems
# for rapidly adding new code) so I'm disabling it this way:
DEVTOOLS_DOCUMENT=devtools::document(roclets=c('namespace', 'rd'))
document: all
	@mkdir -p man
	Rscript -e "library(methods); ${DEVTOOLS_DOCUMENT}"

RSRC = $(shell find ./R)

install: ${PKG}

${PKG}: all ${RSRC}
	R CMD INSTALL --no-test-load .

.PHONY: all doc clean test install
