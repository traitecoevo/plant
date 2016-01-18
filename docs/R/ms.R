## This works around the (frankly) weird misbehaviour of
## rmarkdown::pandoc_convert's wd argument.  I may be misusing it, but
## it seems basically not to work to me.  So instead we'll do the
## chdir ourselves.
pandoc_build <- function(file, template=NULL, format="pdf", ...) {
  if (!grep("\\.md$", file)) {
    stop("Expected a markdown input")
  }

  if (!is.null(template)) {
    template <- normalizePath(template, mustWork=TRUE)
  }
  owd <- setwd(dirname(file))
  on.exit(setwd(owd))

  file_local <- basename(file)
  file_out   <- sub("\\.md$", paste0(".", format), file_local)

  ## Simplify the call a bit, if possible.
  if (is.null(template)) {
    args <- NULL
  } else {
    wd <- getwd()
    if (substr(template, 1, nchar(wd)) == wd) {
      template <- sub("^/+", "",
                      substr(template, nchar(wd) + 1, nchar(template)))
    }
    args <- list(sprintf("--template=%s", template))
  }
  pandoc_convert(file_local, output=file_out, options=args, verbose=FALSE, ...)
}

## TODO: move this into remake
## TODO: options for saying what command is being run
latex_build <- function(filename, bibliography=NULL,
                        chdir=TRUE, interaction="nonstopmode",
                        max_attempts=5L, clean=FALSE, engine="pdflatex") {
  if (chdir && dirname(filename) != "") {
    owd <- setwd(dirname(filename))
    on.exit(setwd(owd))
    filename <- basename(filename)
  }

  res <- run_latex(filename, interaction, engine)
  if(engine=="xelatex") {
      res <- run_latex(filename, interaction, engine)
  }
  if (!is.null(bibliography)) {
    run_bibtex(filename)
    res <- run_latex(filename, interaction, engine)
  }

  pat <- c("Rerun to get cross-references right", # labels
           "Rerun to get citations correct",      # bibtex
           "Rerun to get outlines right")         # tikz
  isin <- function(p, x) {
    any(grepl(p, x))
  }
  for (i in seq_len(max_attempts)) {
    if (any(vapply(pat, isin, logical(1), res))) {
      res <- run_latex(filename, interaction, engine)
    } else {
      break
    }
  }

  if (clean) {
    latex_clean(filename)
  }

  invisible(NULL)
}

latex_clean <- function(filename) {
  filebase <- sub(".tex$", "", filename)
  exts <- c(".log", ".aux", ".bbl", ".blg", ".fls", ".out", ".snm",
            ".nav", ".tdo", ".toc")
  aux <- paste0(filebase, exts)
  file.remove(aux[file.exists(aux)])
}

run_latex <- function(filename, interaction="nonstopmode", engine="pdflatex") {
  args <- c(paste0("-interaction=", interaction),
            "-halt-on-error",
            filename)
  callr::call_system(Sys_which(engine), args)
}

run_bibtex <- function(filename) {
  callr::call_system(Sys_which("bibtex"), sub(".tex$", "", filename))
}

Sys_which <- function(x) {
  ret <- Sys.which(x)
  if (ret == "") {
    stop(sprintf("%s not found in $PATH", x))
  }
  ret
}
