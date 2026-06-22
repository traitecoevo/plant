# new_strategy_scaffolder.R
#
# Scaffold a new plant strategy (and, optionally, environment) by cloning an
# existing model. This generates the boilerplate C++/R files and wires the new
# `<Strategy, Environment>` pair into every place the templated core and the R
# dispatch tables need to learn about it. You still implement the biology by
# hand afterwards (see the `new-strategy` skill, .claude/skills/new-strategy/).
#
# Usage (from the package root):
#
#   source("scripts/new_strategy_scaffolder.R")
#
#   # 1. New strategy WITH its own new environment (clone both from FF16):
#   create_strategy_scaffold("XX99", template_strategy = "FF16")
#
#   # 2. New strategy that REUSES an existing environment (issue #274):
#   #    clone the FF16 biology but share the existing FF16_Environment.
#   #    No `<name>_Environment` files, yml block, or bindings are generated.
#   create_strategy_scaffold("FF16r", template_strategy = "FF16",
#                            environment = "FF16")
#
# After scaffolding: edit the biology, then `make rebuild` and add tests.

library(whisker)
library(yaml)

# ---- small file helpers -----------------------------------------------------

R6_yaml_path <- function() "inst/RcppR6_classes.yml"

updating_message <- function(file) message(sprintf("\t - updating file: %s", file))
creating_message <- function(file) message(sprintf("\t - creating file: %s", file))

read_lines  <- function(file) readLines(file, warn = FALSE)
write_lines <- function(lines, file) writeLines(lines, file)

# Insert `new_lines` immediately after every line containing `anchor`
# (matched as a fixed substring). Errors if the anchor is never found, so a
# stale anchor fails loudly instead of silently producing a broken model.
insert_after <- function(lines, anchor, new_lines, file = "") {
  idx <- which(grepl(anchor, lines, fixed = TRUE))
  if (length(idx) == 0L) {
    stop(sprintf("scaffolder anchor not found in %s:\n  %s", file, anchor))
  }
  for (i in rev(idx)) {
    lines <- append(lines, new_lines, after = i)
  }
  lines
}

# Apply `insert_after` to a file in place.
update_file_insert <- function(file, anchor, new_lines) {
  updating_message(file)
  lines <- read_lines(file)
  lines <- insert_after(lines, anchor, new_lines, file = file)
  write_lines(lines, file)
  invisible(lines)
}

# ---- validation -------------------------------------------------------------

check <- function(name, template_strategy, environment) {
  r6 <- yaml::read_yaml(R6_yaml_path())

  if (length(names(r6)) == 0)
    stop(paste(name, "Cannot find", R6_yaml_path()))
  if (name %in% names(r6))
    stop(paste(name, "is reserved, try again with a different strategy name."))
  if (paste0(name, "_Strategy") %in% names(r6))
    stop(paste("Strategy name:", name, "is already in use, try again with a different strategy name."))
  if (!paste0(template_strategy, "_Strategy") %in% names(r6))
    stop(paste("Template strategy:", template_strategy, "could not be found, please try again."))
  if (!is.null(environment) &&
      !paste0(environment, "_Environment") %in% names(r6))
    stop(paste("Environment to reuse:", environment,
               "could not be found (looked for", paste0(environment, "_Environment"), ")."))
}

# ---- file scaffolding -------------------------------------------------------

# Read a template file, rename the strategy (both `Strategy` and lower-case
# token forms), write to `out_file`. When `reuse_env` is TRUE the environment
# references are rewired back to the reused environment instead of being
# renamed to a (non-existent) `<name>_Environment`.
template_file <- function(name, template_strategy, file, env,
                          out_file = gsub(template_strategy, name,
                                      gsub(tolower(template_strategy), tolower(name), file))) {

  creating_message(out_file)

  comment_char <- ifelse(grepl("R$", file), "#", "//")
  comment <- paste(comment_char, "Built from", file, "on", date(),
                   "using the scaffolder, from the strategy:", template_strategy)

  l_name      <- tolower(name)
  l_strategy  <- tolower(template_strategy)

  out <- read_lines(file)
  out <- gsub(template_strategy, name, out)
  out <- gsub(l_strategy, l_name, out)

  if (env$reuse) {
    # The clone above renamed FF16_Environment -> <name>_Environment etc.
    # Point those references back at the reused environment. Order matters:
    # replace the longer `_Environment` token before the `_Env` token, which
    # is a prefix of it.
    out <- gsub(paste0(name, "_Environment"), env$class, out, fixed = TRUE)
    out <- gsub(paste0(name, "_Env"),         env$type,  out, fixed = TRUE)
    # ... and the header include path.
    out <- gsub(paste0(l_name, "_environment.h"),
                paste0(tolower(env$model), "_environment.h"), out, fixed = TRUE)
  }

  out <- c(comment, out)
  write_lines(out, out_file)
}

scaffold_files <- function(name, template_strategy, env) {
  lt <- tolower(template_strategy)
  files <- c(
    paste0("R/", lt, ".R"),
    paste0("src/", lt, "_strategy.cpp"),
    paste0("src/", lt, "_node.cpp"),
    paste0("inst/include/plant/models/", lt, "_strategy.h"),
    paste0("tests/testthat/test-strategy-", lt, ".R")
  )
  # Only clone the environment header when the new strategy owns its environment.
  if (!env$reuse) {
    files <- c(files,
               paste0("inst/include/plant/models/", lt, "_environment.h"))
  }
  for (file in files) {
    template_file(name, template_strategy, file, env)
  }
}

# ---- RcppR6_classes.yml -----------------------------------------------------

update_classes_yml <- function(name, template_strategy, env) {
  file  <- R6_yaml_path()
  lines <- read_lines(file)

  # 1. Add the new <Strategy, Environment> pair to every `concrete:` block that
  #    currently lists the template pair, plus the IndividualRunner OdeRunner row.
  pair_anchor <- sprintf('["%s": "plant::%s_Strategy", "%s_Env": "plant::%s_Environment"]',
                         template_strategy, template_strategy, template_strategy, template_strategy)
  pair_new <- sprintf('      - ["%s": "plant::%s_Strategy", "%s": "plant::%s"]',
                      name, name, env$type, env$class)
  lines <- insert_after(lines, pair_anchor, pair_new, file = file)

  runner_anchor <- sprintf('["%s": "plant::tools::IndividualRunner<plant::%s_Strategy, plant::%s_Environment>"]',
                           template_strategy, template_strategy, template_strategy)
  runner_new <- sprintf('      - ["%s": "plant::tools::IndividualRunner<plant::%s_Strategy, plant::%s>"]',
                        name, name, env$class)
  lines <- insert_after(lines, runner_anchor, runner_new, file = file)

  # 2. Copy the template's `<Strategy>:` yml block (up to the next blank line)
  #    and append a renamed copy.
  copy_block <- function(lines, block_name, new_block_name) {
    start <- which(grepl(paste0("^", block_name, ":"), lines))[1]
    if (is.na(start)) stop("Could not find yml block: ", block_name)
    rel_blanks <- which(grepl("^\\s*$", lines)) - start
    end <- start + min(rel_blanks[rel_blanks > 0])
    block <- gsub(template_strategy, name, lines[start:(end - 1)])
    block
  }

  updating_message(file)
  appendix <- c(
    "",
    paste("# The following strategy was built from", template_strategy, "on", date()),
    copy_block(lines, paste0(template_strategy, "_Strategy"),
               paste0(name, "_Strategy")))

  # Only add a new environment yml block when the strategy owns its environment.
  if (!env$reuse) {
    appendix <- c(appendix,
      "",
      paste("# The following environment was built from", template_strategy, "on", date()),
      copy_block(lines, paste0(template_strategy, "_Environment"),
                 paste0(name, "_Environment")))
  }

  lines <- c(lines, appendix)
  write_lines(lines, file)
}

# ---- inst/include/plant.h ---------------------------------------------------

update_plant_h <- function(name, template_strategy) {
  lt <- tolower(template_strategy)
  ln <- tolower(name)
  update_file_insert(
    "inst/include/plant.h",
    sprintf("#include <plant/models/%s_strategy.h>", lt),
    sprintf("#include <plant/models/%s_strategy.h>", ln))
}

# ---- src/individual_runner.cpp ----------------------------------------------
# (formerly src/individual_tools.cpp — renamed in the codebase)

update_individual_runner <- function(name, template_strategy, env) {
  file <- "src/individual_runner.cpp"
  updating_message(file)
  block <- whisker.render(
"
// [[Rcpp::export]]
plant::Internals {{name}}_oderunner_individual_internals(
  const odelia::ode::Solver<plant::tools::IndividualRunner<plant::{{name}}_Strategy, plant::{{env_class}}>>& obj) {
  return obj.get_system().individual.r_internals();
}
", list(name = name, env_class = env$class))
  lines <- read_lines(file)
  write_lines(c(lines, strsplit(block, "\n")[[1]]), file)
}

# ---- tests/testthat/helper-plant.R ------------------------------------------

update_test_helper <- function(name, template_strategy, env) {
  file <- "tests/testthat/helper-plant.R"
  updating_message(file)
  lines <- read_lines(file)
  lines <- insert_after(lines,
    sprintf("    %s=%s_Strategy", template_strategy, template_strategy),
    sprintf("    %s=%s_Strategy,", name, name), file = file)
  lines <- insert_after(lines,
    sprintf('    %s="%s_Env"', template_strategy, template_strategy),
    sprintf('    %s="%s",', name, env$type), file = file)
  lines <- insert_after(lines,
    sprintf("    %s=%s_hyperpar", template_strategy, template_strategy),
    sprintf("    %s=%s_hyperpar,", name, name), file = file)
  write_lines(lines, file)
}

# ---- R/strategy_support.R ---------------------------------------------------

update_strategy_support <- function(name, template_strategy, env) {
  file <- "R/strategy_support.R"
  updating_message(file)
  lines <- read_lines(file)

  ins <- function(lines, anchor, new) insert_after(lines, anchor, new, file = file)

  # make_hyperpar()
  lines <- ins(lines, sprintf("         %s=make_%s_hyperpar,", template_strategy, template_strategy),
                      sprintf("         %s=make_%s_hyperpar,", name, name))
  # param_hyperpar()
  lines <- ins(lines, sprintf("         %s_Strategy=%s_hyperpar,", template_strategy, template_strategy),
                      sprintf("         %s_Strategy=%s_hyperpar,", name, name))
  # hyperpar()
  lines <- ins(lines, sprintf("         %s=%s_hyperpar,", template_strategy, template_strategy),
                      sprintf("         %s=%s_hyperpar,", name, name))
  # environment_type()
  lines <- ins(lines, sprintf('         %s=sprintf("%s_Env"),', template_strategy, template_strategy),
                      sprintf('         %s=sprintf("%s"),', name, env$type))
  # Environment(): name -> environment constructor (always).
  lines <- ins(lines, sprintf("         %s=%s_Environment(),", template_strategy, template_strategy),
                      sprintf("         %s=%s(),", name, env$class))
  # Environment(): the "<Env>=" alias only when the strategy owns a NEW env type.
  if (!env$reuse) {
    lines <- ins(lines, sprintf("         %s_Env=%s_Environment(),", template_strategy, template_strategy),
                        sprintf("         %s=%s(),", env$type, env$class))
  }
  # expand_state()
  lines <- ins(lines, sprintf("    %s = %s_expand_state(results),", template_strategy, template_strategy),
                      sprintf("    %s = %s_expand_state(results),", name, name))
  # node_schedule_default()
  lines <- ins(lines,
    sprintf('"Parameters<%s,%s_Env>"=`node_schedule_default__Parameters___%s__%s_Env`',
            template_strategy, template_strategy, template_strategy, template_strategy),
    sprintf('         "Parameters<%s,%s>"=`node_schedule_default__Parameters___%s__%s`,',
            name, env$type, name, env$type))
  # make_node_schedule()
  lines <- ins(lines,
    sprintf('"Parameters<%s,%s_Env>"=`make_node_schedule__Parameters___%s__%s_Env`',
            template_strategy, template_strategy, template_strategy, template_strategy),
    sprintf('         "Parameters<%s,%s>"=`make_node_schedule__Parameters___%s__%s`,',
            name, env$type, name, env$type))

  write_lines(lines, file)
}

# ---- driver -----------------------------------------------------------------

create_strategy_scaffold <- function(name, template_strategy = "FF16",
                                      environment = NULL) {

  reuse <- !is.null(environment)
  env <- list(
    reuse = reuse,
    model = if (reuse) environment else name,
    type  = if (reuse) paste0(environment, "_Env")          else paste0(name, "_Env"),
    class = if (reuse) paste0(environment, "_Environment")  else paste0(name, "_Environment"))

  if (reuse) {
    message(sprintf("Making new strategy '%s' from '%s', reusing '%s' environment",
                    name, template_strategy, env$class))
  } else {
    message(sprintf("Making new strategy '%s' (with its own environment) from '%s'",
                    name, template_strategy))
  }

  check(name, template_strategy, environment)
  scaffold_files(name, template_strategy, env)
  update_classes_yml(name, template_strategy, env)
  update_plant_h(name, template_strategy)
  update_individual_runner(name, template_strategy, env)
  update_test_helper(name, template_strategy, env)
  update_strategy_support(name, template_strategy, env)

  message("\nDone. Next steps:")
  message("  1. Implement the biology in the generated strategy files.")
  message("  2. Run `make rebuild` to regenerate bindings and compile.")
  message("  3. Add tests and run `devtools::test()`.")
  invisible(TRUE)
}
