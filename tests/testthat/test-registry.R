# The model registry (inst/model-registry.yml) is a record, not a check: it
# holds no second implementation of any equation. These tests check that it
# names things that exist and that it leaves nothing out, not that its LaTeX
# agrees with the code.

registry_symbols <- function(symbols) {
  symbols <- unique(stats::na.omit(unlist(symbols)))
  # Template arguments go: the header spells `Species<T,E>::fecundity` and the
  # registry `plant::Species::fecundity`.
  unique(gsub("<[^<>]*>", "", symbols))
}

plant_headers <- function() {
  dir <- system.file("include", "plant", package = "plant")
  files <- list.files(dir, pattern = "[.]h(pp)?$", recursive = TRUE,
                      full.names = TRUE)
  vapply(files, function(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
         character(1))
}

# Whether a qualified name is in the headers, on its class. `plant::Class::m`
# needs `m` in a header that defines Class, or `Class<...>::m` (an out-of-class
# definition) in any header; `plant::f` needs `f` anywhere. A text check, not a
# parse, but one that refuses a member named on the wrong class.
symbol_in_headers <- function(symbol, headers) {
  parts <- strsplit(sub("^plant::", "", symbol), "::", fixed = TRUE)[[1]]
  word <- function(x) sprintf("\\b%s\\b", x)
  if (length(parts) == 1L) {
    return(any(grepl(word(parts), headers, perl = TRUE)))
  }
  cls <- parts[[length(parts) - 1L]]
  member <- parts[[length(parts)]]
  defines <- grepl(sprintf("\\b(class|struct)\\s+%s\\b", cls), headers, perl = TRUE)
  in_class <- any(defines & grepl(word(member), headers, perl = TRUE))
  out_of_class <- any(grepl(sprintf("\\b%s\\s*(<[^;{}()]*>)?\\s*::\\s*~?%s\\b", cls, member),
                            headers, perl = TRUE))
  in_class || out_of_class
}

test_that("the registry reads, and every reference inside it resolves", {
  reg <- model_registry()
  expect_s3_class(reg, "plant_model_registry")
  expect_identical(reg$package, "plant")
  expect_gt(nrow(reg$equations), 0L)
  expect_gt(nrow(reg$slots), 0L)
  # registry_read() refuses an unknown section, class, kind or slot and a read
  # of nothing; this is the positive half: every section is used.
  used <- c(reg$parameters$section, reg$equations$section)
  expect_setequal(names(reg$sections), unique(used))
})

test_that("every C++ symbol the registry names is in plant's headers, on its class", {
  reg <- model_registry()
  headers <- plant_headers()
  symbols <- registry_symbols(list(reg$slots$code, reg$parameters$code,
                                   reg$equations$code, reg$equations$also,
                                   reg$parameters$aliases))
  found <- vapply(symbols, symbol_in_headers, logical(1), headers = headers)
  expect_true(all(found),
              info = paste("not found on its class in inst/include/plant:",
                           toString(symbols[!found])))
  # The check can fail: a real member named on the wrong class is refused.
  expect_false(symbol_in_headers("plant::Patch::survival_individual", headers))
  expect_true(symbol_in_headers("plant::Node::survival_individual", headers))
})

test_that("every field Parameters binds is a parameter or a listed non-parameter", {
  reg <- model_registry()
  bound <- names(scm_base_parameters("FF16"))
  on_parameters <- grepl("^plant::Parameters::", reg$parameters$code)
  field <- ifelse(is.na(reg$parameters$field),
                  sub("^.*::", "", reg$parameters$code), reg$parameters$field)
  covered <- c(field[on_parameters], names(reg$not_parameters))
  expect_setequal(covered, bound)
  expect_false(anyDuplicated(covered) > 0L)
})

test_that("a document takes an equation's LaTeX from the registry", {
  reg <- model_registry()
  id <- reg$equations$id[[1]]
  expect_identical(registry_equation(id, reg),
                   trimws(reg$equations$latex[[1]]))
  expect_error(registry_equation("no_such_equation", reg), "No equation")
})

test_that("the values are read off the live model, not stored", {
  reg <- model_registry()
  v <- registry_values(reg, "FF16")
  p <- scm_base_parameters("FF16")
  i <- match("plant::Parameters::max_patch_lifetime", v$code)
  expect_false(is.na(i))
  expect_identical(v$value[[i]], p$max_patch_lifetime)
})

test_that("the reader refuses what the schema forbids", {
  base <- yaml::read_yaml(system.file("model-registry.yml", package = "plant"))
  refuse <- function(edit, pattern) {
    x <- edit(base)
    f <- withr::local_tempfile(fileext = ".yml")
    yaml::write_yaml(x, f)
    expect_error(model_registry(f), pattern)
  }
  first_eq <- names(base$equations)[[1]]
  refuse(function(x) { x$equations[[first_eq]]$reads <- list("nothing_here"); x },
         "neither a parameter nor a slot")
  refuse(function(x) { x$equations[[first_eq]]$kind <- "guess"; x },
         "unknown kind")
  refuse(function(x) { x$equations[[first_eq]]$section <- "nowhere"; x },
         "unknown section")
  refuse(function(x) { x$slots[[first_eq]] <- x$slots[[1]]; x },
         "ids used twice")
  refuse(function(x) { x$equations[[first_eq]]$latex <- NULL; x },
         "lacks latex")
  refuse(function(x) { x$extra <- 1; x }, "unknown top-level keys")
})
