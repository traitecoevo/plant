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

template <- '
targets:
  maker_scripts.yml:
    depends: "./.scripts"
    command: system("./bootstrap.R")
  scripts:
    depends:
      - maker_scripts.yml
{{#scripts}}
      - {{script}}.md
{{/scripts}}
{{#scripts}}
  {{script}}.Rmd:
    command: sowsear("{{script}}.R", I("Rmd"))
  {{script}}.md:
    knitr: true
{{/scripts}}'

scripts <- sub("\\.R$", "", readLines(".scripts"))
vals <- list(scripts=iteratelist(scripts, value="script"))
yml <- drop_blank(whisker.render(template, vals))
update_file(yml, "maker_scripts.yml")
