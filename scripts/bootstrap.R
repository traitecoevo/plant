#!/usr/bin/Rscript
library(whisker)
update_file <- function(str, filename) {
  if (!file.exists(filename) || !identical(readLines(filename), str)) {
    writeLines(str, filename)
  }
}
drop_blank <- function(x) {
  x <- gsub("\n[[:space:]]*\n", "\n", x)
  gsub("(^\n|\n$)", "", x)
}

load_scripts <- function(filename) {
  scripts <- readLines(filename)
  scripts <- sub("\\s*#.*$", "", scripts)
  scripts <- scripts[!grepl("^\\s*$", scripts)]
  sub("\\.R$", "", scripts)
}

template <- '
targets:
  remake_scripts.yml:
    depends: "./.scripts"
    command: system("./bootstrap.R")
  scripts:
    depends:
      - remake_scripts.yml
{{#scripts}}
      - {{script}}.md
{{/scripts}}
{{#scripts}}
  {{script}}.Rmd:
    command: sowsear("{{script}}.R")
  {{script}}.md:
    knitr: true
{{/scripts}}'

scripts <- load_scripts(".scripts")
vals <- list(scripts=iteratelist(scripts, value="script"))
yml <- drop_blank(whisker.render(template, vals))
update_file(yml, "remake_scripts.yml")
