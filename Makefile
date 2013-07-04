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

document: all
	@mkdir -p man
	Rscript -e "library(methods); devtools::document()"

RSRC = $(shell find ./R)

install: ${PKG}

${PKG}: all ${RSRC}
	R CMD INSTALL --no-test-load .

.PHONY: all doc clean test install
