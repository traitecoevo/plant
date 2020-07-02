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

  # add the extra templates below the FF16 ones
  f <- function(x) {
      switch(x,
      "      - [\"FF16\": \"plant::FF16_Strategy\", \"FF16_Env\": \"plant::FF16_Environment\"]" = c( x, whisker.render(
      "      - [\"{{name}}\": \"plant::{{name}}_Strategy\", \"{{name}}_Env\": \"plant::{{name}}_Environment\"]", list(name=name))),
      "      - [\"FF16\": \"plant::tools::PlantRunner<plant::FF16_Strategy,plant::FF16_Environment>\"]"=c(x,whisker.render(
      "      - [\"{{name}}\": \"plant::tools::PlantRunner<plant::{{name}}_Strategy,plant::{{name}}_Environment>\"]", list(name=name))),
      x)
    }
  update_file(file, name, f) -> r6_templates

  # add the strategy
  str_start <- which(grepl(paste0("^", strategy, "_Strategy"), r6_templates))[1]
  str_candidates <- which(grepl("^$", r6_templates)) # get blank lines
  str_candidates <- str_candidates - str_start
  str_entry <- str_start + min(str_candidates <- str_candidates[str_candidates > 0])

  # grab the parts of the yaml strategy that we will copy and replace
  yml_strategy <- c(
    paste('# The following strategy was built from', strategy, 'on', date()),
    lapply(r6_templates[str_start:str_entry], function (x) {
      gsub(strategy, name, x)
  }))

  # add the environment
  env_start <- which(grepl(paste0("^", strategy, "_Environment"), r6_templates))[1]
  env_candidates <- which(grepl("^$", r6_templates)) # get blank lines
  env_candidates <- env_candidates - env_start
  env_entry <- env_start + min(env_candidates <- env_candidates[env_candidates > 0])
  
  # grab the parts of the yaml strategy that we will copy and replace
  yml_environment <- c(
    paste('# The following environment was built from', strategy, 'on', date()),
    lapply(r6_templates[env_start:env_entry], function (x) {
      gsub(strategy, name, x)
    }))

  append_to_file(file, yml_strategy, str_entry, FALSE)
  append_to_file(file, yml_environment, env_entry, FALSE)
}

# updates src/plant_tools.cpp
update_plant_tools <- function (name) {
  whisker.render("
// [[Rcpp::export]]
plant::Internals {{name}}_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::PlantRunner<plant::{{name}}_Strategy, plant::{{name}}_Environment>>& obj) {
  return obj.obj.plant.r_internals();
}

", list(name=name)) -> debt

  append_to_file("src/plant_tools.cpp", debt, 1)
}

# updates inst/include/plant.h
update_plant <- function (name) {

  # add the extra templates below the FF16 ones
  f <- function(x) {
      if(x != "#include <plant/models/ff16_strategy.h>") return(x)
      c(x, 
        whisker.render("#include <plant/models/{{name}}_strategy.h>",
          list(name=tolower(name))
        ))
    }

  update_file("inst/include/plant.h", name, f)
}

# Updates helper-plant's list of strategies
update_test_helper <- function(name) {
  t1 <- whisker.render("    {{name}}={{name}}_Strategy,", list(name=name))
  t2 <- whisker.render("    {{name}}={{name}}_Environment,", list(name=name))
  t3 <- whisker.render("    {{name}}={{name}}_hyperpar,", list(name=name))
  
  
  # add the extra templates below the FF16 ones
  f <- function(x) {
          switch(x, 
            "    FF16=FF16_Strategy)"=c(t1, x),
            "    FF16=FF16_Environment)"=c(t2, x),
            "    FF16=FF16_hyperpar)"=c(t3, x),
            x)
        }
  update_file("tests/testthat/helper-plant.R", name, f)
}

update_scm_support <- function (name) {
  t1 <- whisker.render("         {{name}}=make_FF16_hyperpar,", list(name=name))
  t2 <- whisker.render("         {{name}}=FF16_hyperpar,", list(name=name))

  # add the extra templates below the FF16 ones
  f <- function(x) {
      switch(x, 
        "         FF16=make_FF16_hyperpar,"=c(t1, x),
        "         FF16=FF16_hyperpar,"=c(t2, x),
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
    paste0('inst/include/plant/models/', tolower(strategy), '_strategy.h'),
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
  update_plant(name)
  update_plant_tools(name)
  update_test_helper(name)
  update_scm_support(name)
}
