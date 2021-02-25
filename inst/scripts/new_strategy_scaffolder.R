library(whisker)
library(yaml)

R6_yaml_path <- function() {
  paste0("inst/RcppR6_classes.yml")
}

check <- function(name, template_strategy) {
  r6 <- yaml::read_yaml(R6_yaml_path())

  if(length(names(r6)) == 0) 
    stop(paste(name, "Cannot find", R6_yaml_path()))
  if(name %in% names(r6)) 
    stop(paste(name, "is reserved, try again with a different strategy name."))
  if(paste0(name, "_Strategy") %in% names(r6))
    stop(paste("Strategy name:", name, "is already in use, try again with a different strategy name."))
  if(! paste0(template_strategy, "_Strategy") %in% names(r6))
    stop(paste("Strategy name:", template_strategy, "could not be found, please try again."))
}

updating_message <- function(file) {
  message(sprintf("\t - updating file: %s", file))
}

creating_message <- function(file) {
  message(sprintf("\t - creating file: %s", file))
}

# update file by applying the function f to each line of the file
# allows for find and replace like modifications
update_file <- function(file, templates, name, template_strategy, sep = "\n") {


  # add the extra templates after the templated ones
  update_txt <- function(template_function) {
    function(x, sep)
      gsub(template_function(template_strategy), 
           sprintf("%s%s%s", 
                   template_function(template_strategy), 
                   sep, 
                   template_function(name)), 
           x, fixed= TRUE)
  }

  updating_message(file)
  readLines(file) -> out

  for(v in templates) {
    template_function <- function(name) {
      whisker.render(v, list(name=name))}
    out <- update_txt(template_function)(out, sep)
  }
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
update_classes_yml <- function (name, template_strategy) {

  file <- R6_yaml_path()

  templates <- c(
    "      - [\"{{name}}\": \"plant::{{name}}_Strategy\", \"{{name}}_Env\": \"plant::{{name}}_Environment\"]",
    "      - [\"{{name}}\": \"plant::tools::IndividualRunner<plant::{{name}}_Strategy, plant::{{name}}_Environment>\"]")

  update_file(file, templates, name, template_strategy) -> r6_templates

  # add the strategy
  str_start <- which(grepl(paste0("^", template_strategy, "_Strategy"), r6_templates))[1]
  str_candidates <- which(grepl("^$", r6_templates)) # get blank lines
  str_candidates <- str_candidates - str_start
  str_entry <- str_start + min(str_candidates <- str_candidates[str_candidates > 0])

  # grab the parts of the yaml strategy that we will copy and replace
  yml_strategy <- c(
    paste('\n# The following strategy was built from', template_strategy, 'on', date()),
    lapply(r6_templates[str_start:str_entry], function (x) {
      gsub(template_strategy, name, x)
  }))

  # add the environment
  env_start <- which(grepl(paste0("^", template_strategy, "_Environment"), r6_templates))[1]
  env_candidates <- which(grepl("^$", r6_templates)) # get blank lines
  env_candidates <- env_candidates - env_start
  env_entry <- env_start + min(env_candidates <- env_candidates[env_candidates > 0])
  
  # grab the parts of the yaml strategy that we will copy and replace
  yml_environment <- c(
    paste('\n# The following environment was built from', template_strategy, 'on', date()),
    lapply(r6_templates[env_start:env_entry], function (x) {
      gsub(template_strategy, name, x)
    }))

  append_to_file(file, c(yml_strategy, yml_environment), env_entry, FALSE)
}

# updates src/plant_tools.cpp
update_plant_tools <- function (name) {
  whisker.render("
// [[Rcpp::export]]
plant::Internals {{name}}_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::IndividualRunner<plant::{{name}}_Strategy, plant::{{name}}_Environment>>& obj) {
  return obj.obj.plant.r_internals();
}

", list(name=name)) -> debt

  append_to_file("src/plant_tools.cpp", debt, 1)
}

# updates inst/include/plant.h
update_plant <- function (name, template_strategy) {

  templates <- c("#include <plant/models/{{name}}_strategy.h>")

  update_file("inst/include/plant.h", templates, tolower(name), tolower(template_strategy))
}

# Updates helper-plant's list of strategies
update_test_helper <- function(name, template_strategy) {

  templates <- c(
  '    {{name}}={{name}}_Strategy',
  '    {{name}}={{name}}_fixed_environment(...)',
  '    {{name}}={{name}}_test_environment(...)',
  '    {{name}}={{name}}_hyperpar'
  )
  
  update_file("tests/testthat/helper-plant.R", templates, name, template_strategy, sep=",\n")
}

update_strategy_support <- function (name, template_strategy) {
  
  templates <- c(
  '         {{name}}=make_{{name}}_hyperpar',
  '         {{name}}={{name}}_hyperpar',
  '         {{name}}_Strategy={{name}}_hyperpar',
  '    {{name}}={{name}}_make_environment(...)',
  '         "Parameters<{{name}},{{name}}_Env>"=`cohort_schedule_max_time_default__Parameters___{{name}}__{{name}}_Env`',
  '         "Parameters<{{name}},{{name}}_Env>"=`cohort_schedule_default__Parameters___{{name}}__{{name}}_Env`',
  '         "Parameters<{{name}},{{name}}_Env>"=`make_cohort_schedule__Parameters___{{name}}__{{name}}_Env`'
  )

  update_file("R/strategy_support.R", templates, name, template_strategy, sep=",\n")
}

# Reads the file, finds and replaces both lower case and uppercase of strategy
# with name and writes the new file to out_file
template_file <- function (name, template_strategy,  file, 
  out_file = gsub(template_strategy, name, 
              gsub(tolower(template_strategy), tolower(name), file))) {

  creating_message(out_file)

  comment_char <- ifelse(grepl("R$", file), '#', '//')
  comment <- paste(comment_char, 'Built from ', file, 'on', date(), 
    'using the scaffolder, from the strategy: ',  template_strategy)

  l_name <- tolower(name)
  l_strategy <- tolower(template_strategy)

  readLines(file) -> raw

  lapply(raw, function(x) gsub(template_strategy, name, x)) -> out
  lapply(out, function(x) gsub(l_strategy, l_name, x)) -> out
  unlist(c(comment, out)) -> out

  writeLines(out, out_file)
} 

scaffold_files <- function (name, template_strategy) {
  files <- c(
    paste0('R/', tolower(template_strategy), '.R'),
    paste0('src/', tolower(template_strategy), '_strategy.cpp'),
    paste0('src/', tolower(template_strategy), '_cohort.cpp'),
    paste0('inst/include/plant/models/', tolower(template_strategy), '_strategy.h'),
    paste0('inst/include/plant/models/', tolower(template_strategy), '_environment.h'),
    paste0('tests/testthat/test-strategy-', tolower(template_strategy), '.R')
  )
  for (file in files) {
    template_file(name, template_strategy, file)
  }  
}

create_strategy_scaffold <- function(name, template_strategy="FF16") {
  message(sprintf("Making new strategy with name %s from %s", name, template_strategy))
  check(name, template_strategy)
  scaffold_files(name, template_strategy)
  update_classes_yml(name, template_strategy)
  update_plant(name, template_strategy)
  update_plant_tools(name)
  update_test_helper(name, template_strategy)
  update_strategy_support(name, template_strategy)
}
