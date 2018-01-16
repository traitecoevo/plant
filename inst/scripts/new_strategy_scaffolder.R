library(whisker)
library(yaml)

R6_yaml_path <- function() {
  paste0("inst/RcppR6_classes.yml")
}

check <- function(name, strategy) {
  r6 <- yaml::read_yaml(R6_yaml_path())

  if(length(names(r6)) == 0) 
    stop(paste(name, "Cannot find", R6_yaml_path()))
  if(name %in% names(r6)) 
    stop(paste(name, "is reserved, try again with a different strategy name."))
  if(paste0(name, "_Strategy") %in% names(r6))
    stop(paste("Strategy name:", name, "is already in use, try again with a different strategy name."))
  if(! paste0(strategy, "_Strategy") %in% names(r6))
    stop(paste("Strategy name:", strategy, "could not be found, please try again."))
}

updating_message <- function(file) {
  message(sprintf("\t - updating file: %s", file))
}

creating_message <- function(file) {
  message(sprintf("\t - creating file: %s", file))
}

# update file by applying the function f to each line of the file
# allows for find and replace like modifications
update_file <- function(file, name, f) {

  updating_message(file)
  readLines(file) -> raw
  lapply(raw, f) -> out
  paste0(unlist(out)) -> out  
  writeLines(out, file)
  invisible(out)
}

# append text to end of file
append_to_file <- function(file, debt, after, verbose=TRUE) {

  if(verbose)
    updating_message(file)
  readLines(file) -> raw
  append(raw, unlist(debt), length(raw) + after) -> out
  writeLines(out, file)
}

# make the RccpR6 file
update_classes_yml <- function (name, strategy) {

  file <- R6_yaml_path()

  # add the extra templates below the FF16r ones
  f <- function(x) {
      if(x != "      - [\"FF16r\": \"plant::FF16r_Strategy\"]") return(x)
      c(x, 
        whisker.render(
          "      - [\"{{name}}\": \"plant::{{name}}_Strategy\"]",
          list(name=name)
        ))
    }
  update_file(file, name, f) -> r6_templates

  # add the strategy

  start_point <- which(grepl(paste0("^", strategy, "_Strategy"), r6_templates))[1]
  candidates <- which(grepl("^$", r6_templates)) # get blank lines
  candidates <- candidates - start_point
  entry <- start_point +  min(candidates <- candidates[ candidates > 0])

  # grab the parts of the yaml strategy that we will copy and replace
  yml_strategy <- c(paste('# The following strategy was built from', strategy, 'on', date()),
  lapply(r6_templates[start_point:entry], function (x) {
    gsub(strategy, name, x)
  }))

  append_to_file(file, yml_strategy, entry, FALSE)
}

# updates src/plant_tools.cpp
update_plant_tools <- function (name) {
  whisker.render("
// Template found in scripts/strategy/new_strategy.R
// [[Rcpp::export]]
double {{name}}_lcp_whole_plant(plant::PlantPlus<plant::{{name}}_Strategy> p) {
  return plant::tools::lcp_whole_plant(p);
}
", list(name=name)) -> debt

  append_to_file("src/plant_tools.cpp", debt, 1)
}

# updates src/plant_plus.cpp
update_plant_plus <- function (name) {

  f <- function(x) {
      if(x != "#include <plant/ff16r_strategy.h>") return(x)
      c(x, 
        whisker.render("#include <plant/{{name}}_strategy.h>",
          list(name=tolower(name))
        ))
    }
  update_file("src/plant_plus.cpp", name, f) 

  # add the technical debt to the end of the file
  whisker.render("
// Template found in scripts/strategy/new_strategy.R
// [[Rcpp::export]]
plant::PlantPlus<plant::{{name}}_Strategy>
{{name}}_plant_to_plant_plus(plant::Plant<plant::{{name}}_Strategy> p,
                          SEXP environment) {
  return plant::plant_to_plant_plus(p, environment);
}
", list(name=name)) -> debt

  append_to_file("src/plant_plus.cpp", debt, 1, FALSE)
}

# updates inst/include/plant.h
update_plant <- function (name) {

  # add the extra templates below the FF16r ones
  f <- function(x) {
      if(x != "#include <plant/ff16r_strategy.h>") return(x)
      c(x, 
        whisker.render("#include <plant/{{name}}_strategy.h>",
          list(name=tolower(name))
        ))
    }

  update_file("inst/include/plant.h", name, f)
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
  append_to_file("R/plant.R", list(t1, t2), 5)
}

# Updates helper-plant's list of strategies
update_test_helper <- function(name) {
  t1 <- whisker.render("       {{name}}={{name}}_Strategy,", list(name=name))
  t2 <- whisker.render("       {{name}}={{name}}_hyperpar,", list(name=name))

  # add the extra templates below the FF16r ones
  f <- function(x) {
          switch(x, 
            "       FF16r=FF16r_Strategy)"=c(t1, x),
            "       FF16r=FF16r_hyperpar)"=c(t2, x),
            x)
        }
  update_file("tests/testthat/helper-plant.R", name, f)
}

update_scm_support <- function (name) {
  t1 <- whisker.render("         {{name}}=make_FF16_hyperpar,", list(name=name))
  t2 <- whisker.render("         {{name}}=FF16_hyperpar,", list(name=name))

  # add the extra templates below the FF16r ones
  f <- function(x) {
      switch(x, 
        "         FF16r=make_FF16_hyperpar,"=c(t1, x),
        "         FF16r=FF16_hyperpar,"=c(t2, x),
        x
      )
    }
  update_file("R/scm_support.R", name, f)
}

# Reads the file, finds and replaces both lower case and uppercase of strategy
# with new_name and writes the new file to out_file
template_file <- function (new_name, strategy, file, 
  out_file = gsub(strategy, new_name, 
              gsub(tolower(strategy), tolower(new_name), file)),
  test = TRUE) {

  ar <- function (..., test = FALSE) paste0(ifelse(test, './results/',  "./"), ...)
  
  if (test) out_file <- gsub('/', '__', out_file)

  creating_message(out_file)

  comment_char <- ifelse(grepl("R$", file), '#', '//')
  comment <- paste(comment_char, 'Built from ', file, 'on', date(), 
    'using the scaffolder, from the strategy: ',  strategy)

  l_new_name <- tolower(new_name)
  l_strategy <- tolower(strategy)

  readLines(ar(file)) -> raw

  lapply(raw, function(x) gsub(strategy, new_name, x)) -> out
  lapply(out, function(x) gsub(l_strategy, l_new_name, x)) -> out
  unlist(c(comment, out)) -> out

  writeLines(out, ar(out_file, test = test))
} 

scaffold_files <- function (new_name, strategy, test = FALSE) {
  files <- c(
    paste0('R/', tolower(strategy), '.R'),
    paste0('src/', tolower(strategy), '_strategy.cpp'),
    paste0('inst/include/plant/', tolower(strategy), '_strategy.h'),
    paste0('tests/testthat/test-strategy-', tolower(strategy), '.R')
  )
  for (file in files) {
    template_file(new_name, strategy, file, test = test)
  }  
}

create_strategy_scaffold <- function(name, template="FF16") {
  message(sprintf("Making new strategy with name %s from %s", name, template))
  check(name, template)
  scaffold_files(name, template)
  update_classes_yml(name, template)
  update_plant_plus(name)
  update_plant(name)
  update_plant_tools(name)
  update_plant_r(name)
  update_test_helper(name)
  update_scm_support(name)
}
