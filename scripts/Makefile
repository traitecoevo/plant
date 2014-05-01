RSCRIPT_PKGS := $(shell Rscript -e 'library(methods);writeLines(Sys.getenv("R_DEFAULT_PACKAGES"))')
RSCRIPT = Rscript --default-packages="${RSCRIPT_PKGS},methods"

SCRIPTS = $(shell cat .wiki_scripts)
TARGETS = $(patsubst %.R,%.md, ${SCRIPTS})
TARGETS_RMD = $(patsubst %.R,%.Rmd, ${SCRIPTS})

all: ${TARGETS}

%.Rmd: %.R
	${RSCRIPT} -e "library(sowsear); sowsear('$<', 'Rmd')"
%.md: %.Rmd
	./.knit.sh $<
%.pdf: %.md
	pandoc $< -o $@
%.docx: %.md
	pandoc $< -o $@
%.html: %.md
	${RSCRIPT} -e "library(markdown);\
	  opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images');\
	markdownToHTML('$<', '$@', options=opts)"

clean:
	rm -f ${TARGETS} ${TARGETS_RMD}
	rm -rf figure cache

WIKI_DIR = ../wiki
wiki_update:
	./.wiki_update.sh
wiki_publish: wiki_update
	git --git-dir=${WIKI_DIR}/.git push
wiki_reset:
	git --git-dir=${WIKI_DIR}/.git --work-tree=${WIKI_DIR} reset --hard

.PHONY: all clean wiki_update wiki_publish wiki_reset
.SECONDARY:
