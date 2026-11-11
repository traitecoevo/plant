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
  shared <- names(reg$not_parameters)[!grepl(".", names(reg$not_parameters), fixed = TRUE)]
  covered <- c(field[on_parameters], shared)
  expect_setequal(covered, bound)
  expect_false(anyDuplicated(covered) > 0L)
})

test_that("every field of each model's pars is a parameter or a listed non-parameter", {
  reg <- model_registry()
  for (m in reg$models$id) {
    bound <- names(get(paste0(m, "_Model"))()$pars)
    lineage <- plant:::registry_lineage(reg$models, m)
    own <- reg$parameters[!is.na(reg$parameters$model) &
                            reg$parameters$model %in% lineage, ]
    field <- ifelse(is.na(own$field), sub("^.*::", "", own$code), own$field)
    np <- names(reg$not_parameters)
    listed <- unlist(lapply(lineage, function(l)
      sub(paste0("^", l, "[.]"), "", np[startsWith(np, paste0(l, "."))])))
    covered <- c(field, listed)
    expect_setequal(covered, bound)
    expect_false(anyDuplicated(covered) > 0L, info = m)
  }
})

test_that("every model fills every per-model slot, or says why not", {
  reg <- model_registry()
  per_model <- reg$slots$id[!(reg$slots$per_model %in% "FALSE")]
  entries <- rbind(reg$parameters[c("id", "model", "fills")],
                   reg$equations[c("id", "model", "fills")])
  for (i in seq_len(nrow(reg$models))) {
    m <- reg$models$id[[i]]
    lineage <- plant:::registry_lineage(reg$models, m)
    filled <- entries$fills[!is.na(entries$model) & entries$model %in% lineage]
    unfilled <- names(reg$models$unfilled[[i]])
    missing <- setdiff(per_model, c(filled, unfilled))
    expect_true(length(missing) == 0L,
                info = sprintf("%s fills no %s", m, toString(missing)))
    expect_length(intersect(filled, unfilled), 0L)
    expect_false(anyDuplicated(filled) > 0L,
                 info = sprintf("%s fills %s twice", m, toString(filled[duplicated(filled)])))
  }
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
  # A model's own parameter is read off its default strategy.
  j <- match("FF16.lma", v$id)
  expect_false(is.na(j))
  expect_identical(v$value[[j]], FF16_Model()$pars$lma)
  # and another model's is not this model's to give.
  k <- match("K93.b_0", v$id)
  expect_false(is.na(k))
  expect_true(is.na(v$value[[k]]))
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
  first_par <- names(base$parameters)[[1]]
  refuse(function(x) { x$parameters[["FF16.made_up"]] <- x$parameters[[first_par]]; x },
         "named like a model's entry but has no model")
  refuse(function(x) {
    e <- x$parameters[[first_par]]; e$model <- "FF16"
    x$parameters[["K93.made_up"]] <- e; x
  }, "so its id starts 'FF16.'")
  refuse(function(x) {
    e <- x$parameters[[first_par]]; e$model <- "NOPE"
    x$parameters[["NOPE.made_up"]] <- e; x
  }, "unknown model 'NOPE'")
  refuse(function(x) { x$models$K93$extends <- "NOPE"; x }, "extends unknown model")
  refuse(function(x) { x$models$K93$extends <- "FF16"; x$models$FF16$extends <- "K93"; x },
         "extends itself")
  refuse(function(x) { x$parameters[[first_par]]$fills <- "nonsense"; x },
         "fills unknown slot")
  refuse(function(x) { x$parameters[[first_par]]$fills <- "g"; x },
         "is shared but fills per-model slot")
  refuse(function(x) { x$parameters[["K93.eta"]]$fills <- "Q"; x },
         "which every model shares")
  refuse(function(x) { x$not_parameters[["FF61.foo"]] <- "typo"; x },
         "names an unknown model")
  refuse(function(x) { x$models$K93$unfilled$w <- ""; x }, "with no reason")
})
