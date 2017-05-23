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



##' Function imported from callr package;  makes it easy to call a
##' system command from R and have it behave.
##'
##' This function uses \code{system2} to call a system command fairly
##' portably.  What it adds is a particular way of dealing with
##' errors.  \code{call_system} runs the command \code{command} with
##' arguments \code{args} (and with optionally set environment
##' variables \code{env}) and hides \emph{all} produced output to
##' stdout and stderr.  If the command fails (currently any nonzero
##' exit code is counted as a failure) then \code{call_system} will
##' throw an R error giving
##' \itemize{
##' \item the full string of the command run
##' \item the exit code of the command
##' \item any \code{errmsg} attribute that might have been returned
##' \item all output that the program produced to either stdout and
##' stderr
##' }
##'
##' This means that a successful invocation of a program produces no
##' output while the unsuccessful invocation throws an error and
##' prints all information to the screen (though this is delayed until
##' failure happens).
##'
##'
##' \code{call_system} also returns the contents of both stderr and
##' stdout \emph{invisibly} so that it can be inspected if needed.
##'
##' The function \code{run_system} does the same thing and will be
##' removed as soon as code that depends on it is out of use.
##'
##' @title Run a system command, stopping on error
##' @param command The system command to be invoked, as a character
##' string.  \code{\link{Sys.which}} is useful here.
##' @param args A character vector of arguments to \code{command}
##' @param env A character vector of name=value pairs to be set as
##' environment variables (see \code{\link{system2}}).
##' @param max_lines Maximum number of lines of program output to
##' print with the error message.  We may prune further to get the
##' error message under \code{getOption("warn.length")}, however.
##' @param p Fraction of the error message to show from the tail of
##' the output if truncating on error (default is 20\% lines are head,
##' 80\% is tail).
##' @param stdout,stderr Passed to \code{system2}.  Set one of these
##' to \code{FALSE} to avoid capturing output from that stream.  Setting
##' both to \code{FALSE} is not recommended.
##' @export
##' @author Rich FitzJohn
call_system <- function(command, args, env=character(), max_lines=20,
                        p=0.8, stdout=TRUE, stderr=TRUE) {
  res <- suppressWarnings(system2(command, args,
                                  env=env, stdout=stdout, stderr=stderr))
  ok <- attr(res, "status")
  if (!is.null(ok) && ok != 0) {
    max_nc <- getOption("warning.length")

    cmd <- paste(c(env, shQuote(command), args), collapse = " ")
    msg <- sprintf("Running command:\n  %s\nhad status %d", cmd, ok)
    errmsg <- attr(cmd, "errmsg")
    if (!is.null(errmsg)) {
      msg <- c(msg, sprintf("%s\nerrmsg: %s", errmsg))
    }
    sep <- paste(rep("-", getOption("width")), collapse="")

    ## Truncate message:
    if (length(res) > max_lines) {
      n <- ceiling(max_lines * p)
      res <- c(head(res, ceiling(max_lines - n)),
               sprintf("[[... %d lines dropped ...]]", length(res) - max_lines),
               tail(res, ceiling(n)))
    }

    ## compute the number of characters so far, including three new lines:
    nc <- (nchar(msg) + nchar(sep) * 2) + 3
    i <- max(1, which(cumsum(rev(nchar(res) + 1L)) < (max_nc - nc)))
    res <- res[(length(res) - i + 1L):length(res)]
    msg <- c(msg, "Program output:", sep, res, sep)
    stop(paste(msg, collapse="\n"))
  }
  invisible(res)
}

