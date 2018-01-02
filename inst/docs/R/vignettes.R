patch_md <- function(filename, dest) {
  'title: "%s"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %%\\VignetteIndexEntry{%s}
  %%\\VignetteEngine{knitr::rmarkdown}
  \\usepackage[utf8]{inputenc}' -> template
  dat <- readLines(filename)
  i <- grep("^---$", dat)
  re <- '^title:\\s+"(.*)"\\s*'
  if (!identical(i, c(1L, 3L)) || !grepl(re, dat[[2]])) {
    stop("malformed markdown file")
  }
  title <- sub(re, "\\1", dat[[2]])
  dat[[2]] <- sprintf(template, title, title)
  writeLines(dat, dest)
}

build_vignettes <- function(cleanup=FALSE) {
  oo <- options(warnPartialMatchArgs=FALSE)
  if (!is.null(oo$warnPartialMatchArgs)) {
    options(oo)
  }
  dest <- "vignettes"
  copy <- c("figures", "figure",
            "mee.bst", "suppmat.sty", "refs.bib",
            "growth_model_pars_core.tex", "growth_model_pars_hyper.tex")
  unlink(file.path(dest, copy), recursive=TRUE)
  if (!cleanup) {
    file.copy(copy, dest, recursive=TRUE, overwrite=TRUE)
    devtools::build_vignettes()
  }
}

cleanup_vignettes <- function() {
  build_vignettes(TRUE)
}
