#!/bin/sh

# Compile Rmd -> md, setting some options.  The options that I want
# are:
#
#   * cache and figure directory names based on the current script so
#     that we don't get collisions between chunks from different files
#   * halt on errors, rather than just keeping on going.
#   * don't reformat code

# First, work out the prefix by stripping the extension.
TARGET=$1

# This check might be worthwhile, because otherwise we generate some
# funky filenames.
if ! echo $TARGET | grep --quiet '.Rmd$'
then
    echo "Error: Does not look like an Rmd file, aborting"
    exit 1
fi

# The basename here will hopefully mean that knitr will do the right
# thing when run in a directory other than that containing the file,
# but this is untested.
PREFIX=$(echo `basename $TARGET` | sed 's/.Rmd$//')

# Common options that I want for everything (must be at least one
# option, but all but fig.height I *always* want.
KNITR_OPTS="error=FALSE, tidy=FALSE, fig.height=5"
KNITR_PATHS="fig.path='figure/${PREFIX}_', cache.path='cache/${PREFIX}_'"
KNITR_OPTS_CMD="knitr::opts_chunk\$set(${KNITR_OPTS}, ${KNITR_PATHS})"

Rscript -e "library(methods); ${KNITR_OPTS_CMD}; knitr::knit('${TARGET}')"
