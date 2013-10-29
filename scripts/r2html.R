#!/usr/bin/env Rscript
file <- commandArgs(TRUE)[[1]]
if (!file.exists(file))
  stop("Cannot see file ", dQuote(file))

if (!require(sowsear))
  stop("Install sowsear from https://github.com/richfitz/sowsear")

file.base <- tools::file_path_sans_ext(file)
file.Rmd  <- paste0(file.base, ".Rmd")
file.md   <- paste0(file.base, ".md")
fig.path  <- paste0("_figures/", file.base, "_")

opts_chunk$set(tidy=FALSE, fig.height=5, fig.path=fig.path)

sowsear(file, "Rmd")  # Generate Rmd file
knit2html(file.Rmd)   # Convert Rmd -> md and md -> html
file.remove(c(file.Rmd, file.md)) # Delete temporary files
unlink(dirname(fig.path), recursive=TRUE)
unlink("cache", recursive=TRUE)
