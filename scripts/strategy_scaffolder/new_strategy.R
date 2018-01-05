#!/bin/usr/RScript

library(whisker)
library(yaml)

root <-  "./../.."
R6_yaml_path <- paste0(root, "/inst/RcppR6_classes.yml")

check <- function(name) {
  r6 <- yaml::read_yaml(R6_yaml_path)

  if(length(names(r6)) == 0) 
    stop(paste(name, "Cannot find", R6_yaml_path))
  if(name %in% names(r6)) 
    stop(paste(name, "is reserved, try again with a different strategy name."))
  if(paste0(name, "_Strategy") %in% names(r6))
    stop(paste("Strategy name:", name, "is allready in use, try again with a different strategy name."))
}


read_template <- function (filename) {
  readChar(filename, file.info(filename)$size)
}

# Makes the R file:
# R/NAME.R
render_r <- function(name) {
  template <- read_template("./templates/NAME.Rtemplate")
  writeLines(
    whisker.render(template, list(name=name)),
    paste0(root, "/R/", tolower(name), ".R")
  )
}

# Makes the new file
# inst/include/plant/NAME_strategy.h
render_h <- function (name) {
  template <- read_template("./templates/NAME_strategy.htempate")
  writeLines(
    whisker.render(template, list(name=name)),
    paste0(root, "/inst/include/plant/", tolower(name), "_strategy.h")
  )
}
# * src/NAME_strategy.cpp
render_cpp <- function (name) {
  template <- read_template("./templates/NAME_strategy.cpptemplate")
  writeLines(
    whisker.render(template, list(name=name,lname=tolower(name))),
    paste0(root, "/src/", tolower(name), "_strategy.cpp")
  )
} 

# make the RccpR6 file
update_classes_yml <- function (name) {
  readLines(R6_yaml_path) -> raw
  # add the extra templates below the FF16r ones
  lapply(raw, function(x) {
      if(x != "      - [\"FF16r\": \"plant::FF16r_Strategy\"]") return(x)
      c(x, 
        whisker.render(
          "      - [\"{{name}}\": \"plant::{{name}}_Strategy\"]",
          list(name=name)
        ))
    }) -> r6_templates
  
  unlist(r6_templates) -> r6_templates

  # add the strategy

  start_point <- which(grepl("^FF16r_Strategy", r6_templates))[1]
  candidates <- which(grepl("^$", r6_templates)) # get blank lines
  candidates <- candidates - start_point
  entry <- start_point +  min(candidates <- candidates[ candidates > 0])

  read_template('./templates/R6Strategy.ymltemplate') -> yml_template
  whisker.render(yml_template, list(name=name)) -> yml_strategy
  append(r6_templates, yml_strategy, entry) -> out

  # write
  paste0(unlist(out)) -> out
  writeLines(
    out, 
    paste0(root, "/inst/RcppR6_classes.yml")
  )
}

# updates src/plant_tools.cpp
update_plant_tools <- function (name) {
  readLines(paste0(root, "/src/plant_tools.cpp")) -> raw
    whisker.render("
// Template found in scripts/strategy/new_strategy.R
// [[Rcpp::export]]
double {{name}}_lcp_whole_plant(plant::PlantPlus<plant::{{name}}_Strategy> p) {
  return plant::tools::lcp_whole_plant(p);
}
", list(name=name)) -> debt

  append(raw, debt, length(raw) + 1) -> out
  # write
  paste0(unlist(out)) -> out
  writeLines(
    out,
    paste0(root, "/src/plant_tools.cpp")
  )
}

# updates src/plant_plus.cpp
update_plant_plus <- function (name) {

  readLines(paste0(root, "/src/plant_plus.cpp")) -> raw
  # add the extra templates below the FF16r ones
  lapply(raw, function(x) {
      if(x != "#include <plant/ff16r_strategy.h>") return(x)
      c(x, 
        whisker.render("#include <plant/{{name}}_strategy.h>",
          list(name=tolower(name))
        ))
    }) -> with_includes

  unlist(with_includes) -> with_includes

  # add the the technical debt to the end of the file

  whisker.render("
// Template found in scripts/strategy/new_strategy.R
// [[Rcpp::export]]
plant::PlantPlus<plant::{{name}}_Strategy>
{{name}}_plant_to_plant_plus(plant::Plant<plant::{{name}}_Strategy> p,
                          SEXP environment) {
  return plant::plant_to_plant_plus(p, environment);
}
", list(name=name)) -> debt

  append(with_includes, debt, length(raw) + 1) -> out
  # write
  paste0(unlist(out)) -> out
  writeLines(
    out, 
    paste0(root, "/src/plant_plus.cpp")
  )
}

# updates inst/include/plant.h
update_plant <- function (name) {
    readLines(paste0(root, "/inst/include/plant.h")) -> raw
  # add the extra templates below the FF16r ones
  lapply(raw, function(x) {
      if(x != "#include <plant/ff16r_strategy.h>") return(x)
      c(x, 
        whisker.render("#include <plant/{{name}}_strategy.h>",
          list(name=tolower(name))
        ))
    }) -> out
  paste0(unlist(out)) -> out
  writeLines(
    out, 
    paste0(root, "/inst/include/plant.h")
  )
}

update_plant_r <- function (name) {
  t1 <- whisker.render("# The following two functions were added by the scaffolder in
# /scripts/strategy_scaffolder/new_strategy.R
##' @export
`plant_to_plant_plus.Plant<{{name}}>` <- function(x, ...) {
  {{name}}_plant_to_plant_plus(x, ...)
}", list(name=name))
  t2 <- whisker.render("##' @export
`lcp_whole_plant.PlantPlus<{{name}}>` <- function(p, ...) {
  {{name}}_lcp_whole_plant(p, ...)
}", list(name=name))

  # add both t1 and t2 at the end of the file

  readLines(paste0(root, "/R/plant.R")) -> raw
  # add the extra templates below the FF16r ones
  append(raw, list(t1, t2), length(raw) + 5) -> out
  paste0(unlist(out)) -> out
  writeLines(
    out, 
    paste0(root, "/R/plant.R")
  )
}

# Builds the tests to tests/testthat/test-strategy-name.R
render_tests <- function(name) {
  template <- read_template("./templates/test-strategy-NAME.Rtemplate")
  writeLines(
    whisker.render(template, list(name=name)),
    paste0(root, "/tests/testthat/test-strategy-", tolower(name), ".R")
  )
}
# Updates helper-plant's list of strategies
update_test_helper <- function(name) {
  t1 <- whisker.render("       {{name}}={{name}}_Strategy,", list(name=name))
  t2 <- whisker.render("       {{name}}={{name}}_hyperpar,", list(name=name))

  readLines(paste0(root, "/tests/testthat/helper-plant.R")) -> raw
  # add the extra templates below the FF16r ones
  lapply(raw, function(x) {
      switch(x, 
        "       FF16r=FF16r_Strategy)"=c(t1, x),
        "       FF16r=FF16r_hyperpar)"=c(t2, x),
        x
      )
    }) -> out
  paste0(unlist(out)) -> out  
  writeLines(out,
    paste0(root, "/tests/testthat/helper-plant.R")
  )
}

update_scm_support <- function (name) {
  t1 <- whisker.render("         {{name}}=make_FF16_hyperpar,", list(name=name))
  t2 <- whisker.render("         {{name}}=FF16_hyperpar,", list(name=name))

  readLines(paste0(root, "/R/scm_support.R")) -> raw
  # add the extra templates below the FF16r ones
  lapply(raw, function(x) {
      switch(x, 
        "         FF16r=make_FF16_hyperpar,"=c(t1, x),
        "         FF16r=FF16_hyperpar,"=c(t2, x),
        x
      )
    }) -> out
  paste0(unlist(out)) -> out
  writeLines(out,
    paste0(root, "/R/scm_support.R")
  )
}
new_strategy <- function(name) {
  check(name)
  render_r(name)
  render_h(name)
  render_cpp(name)
  update_classes_yml(name)
  update_plant_plus(name)
  update_plant(name)
  update_plant_tools(name)
  update_plant_r(name)
  render_tests(name)
  update_test_helper(name)
  update_scm_support(name)
}