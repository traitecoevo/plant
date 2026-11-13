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
  bound <- names(base_parameters("FF16"))
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
    # Only what the model's _Pars holds: a model's constants on its class and
    # its environment's parameters are entries too, but not pars fields.
    own <- reg$parameters[!is.na(reg$parameters$model) &
                            reg$parameters$model %in% lineage &
                            grepl("_Pars::", reg$parameters$code), ]
    field <- ifelse(is.na(own$field), sub("^.*::", "", own$code), own$field)
    np <- names(reg$not_parameters)
    listed <- unlist(lapply(lineage, function(l)
      sub(paste0("^", l, "[.]"), "", np[startsWith(np, paste0(l, "."))])))
    covered <- c(field, listed)
    expect_setequal(covered, bound)
    expect_false(anyDuplicated(covered) > 0L, info = m)
  }
})

test_that("every field a subclass adds beside pars is a parameter or a listed non-parameter", {
  # A subclass's own traits (TF24f's k_acclim, TF24t's thermal set) sit at the
  # strategy's top level, beside the parent's pars, where the pars test cannot
  # see them.
  reg <- model_registry()
  sub <- reg$models[!is.na(reg$models$extends), ]
  expect_gt(nrow(sub), 0L)
  for (i in seq_len(nrow(sub))) {
    m <- sub$id[[i]]
    s <- get(paste0(m, "_Model"))()
    added <- setdiff(names(s), names(get(paste0(sub$extends[[i]], "_Model"))()))
    added <- added[vapply(added, function(f) is.atomic(s[[f]]), logical(1))]
    own <- reg$parameters[reg$parameters$model %in% m &
                            grepl(paste0("^plant::", m, "_Model::"), reg$parameters$code), ]
    field <- ifelse(is.na(own$field), sub("^.*::", "", own$code), own$field)
    np <- names(reg$not_parameters)
    listed <- sub(paste0("^", m, "[.]"), "", np[startsWith(np, paste0(m, "."))])
    missing <- setdiff(added, c(field, listed))
    expect_true(length(missing) == 0L,
                info = sprintf("%s adds %s, which the registry does not cover", m, toString(missing)))
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
    filled <- filled[!is.na(filled)]
    unfilled <- names(reg$models$unfilled[[i]])
    missing <- setdiff(per_model, c(filled, unfilled))
    expect_true(length(missing) == 0L,
                info = sprintf("%s fills no %s", m, toString(missing)))
    expect_length(intersect(filled, unfilled), 0L)
    expect_false(anyDuplicated(filled) > 0L,
                 info = sprintf("%s fills %s twice", m, toString(filled[duplicated(filled)])))
  }
})

test_that("every link resolves in the phylloptim plant is built against", {
  skip_if_not_installed("phylloptim")
  # The pin is checked first: a link is checked against the installed registry,
  # which is the right one only when it is the version plant compiles.
  pinned <- plant:::registry_pinned_version("phylloptim")
  expect_false(is.na(pinned))
  expect_identical(as.character(utils::packageVersion("phylloptim")), pinned)

  reg <- model_registry()
  expect_gt(nrow(reg$links), 0L)
  expect_true(all(reg$links$resolved))
  expect_false(anyNA(reg$links$kind))
  # A linked parameter takes its description from its target, not from here.
  target <- yaml::read_yaml(system.file("model-registry.yml", package = "phylloptim"))
  i <- match("TF24.stem_P50", reg$parameters$id)
  expect_identical(reg$parameters$link[[i]], "phylloptim:stem_P50")
  expect_identical(reg$parameters$symbol[[i]], target$parameters$stem_P50$symbol)
  expect_identical(reg$parameters$units[[i]], target$parameters$stem_P50$units)
  # Every linked parameter still names the plant member that hands it over.
  linked <- reg$parameters[!is.na(reg$parameters$link), ]
  expect_false(anyNA(linked$symbol))
  expect_true(all(startsWith(linked$code, "plant::")))
  # and a page can render a linked equation through the same call.
  expect_identical(registry_equation("phylloptim:carbon", reg),
                   trimws(target$equations$carbon$latex))
})

test_that("the reader refuses a link it cannot resolve", {
  skip_if_not_installed("phylloptim")
  base <- yaml::read_yaml(system.file("model-registry.yml", package = "plant"))
  refuse <- function(edit, pattern) {
    f <- withr::local_tempfile(fileext = ".yml")
    yaml::write_yaml(edit(base), f)
    expect_error(suppressWarnings(model_registry(f)), pattern)
  }
  refuse(function(x) { x$parameters$TF24.stem_P50$link <- "phylloptim:no_such_id"; x },
         "does not hold")
  refuse(function(x) { x$parameters$TF24.stem_P50$link <- "phylloptim:carbon"; x },
         "is an equation in phylloptim's registry, not a parameter")
  refuse(function(x) { x$equations$TF24.assimilation$reads <- list("phylloptim:no_such_id"); x },
         "does not hold")
  # The equivalent mutation of a live entry, so the mechanism is seen to bite.
  refuse(function(x) { x$equations$TF24.assimilation$reads <- list("phylloptim:carbon_kept"); x },
         "does not hold")
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
  p <- base_parameters("FF16")
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
  expect_true(is.na(v$value[[k]]))  # A subclass's trait beside pars is read off the strategy itself, and the
  # parent's parameters off its pars.
  vt <- registry_values(reg, "TF24t")
  expect_identical(vt$value[[match("TF24t.k_i", vt$id)]], TF24t_Model()$k_i)
  expect_identical(vt$value[[match("TF24.lma", vt$id)]], TF24t_Model()$pars$lma)
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
  # Links: to a declared package only, and a linked parameter restates nothing.
  refuse(function(x) { x$parameters$TF24.stem_P50$link <- "leafy:stem_P50"; x },
         "for a package under links")
  refuse(function(x) { x$parameters$TF24.stem_P50$symbol <- "P"; x },
         "must not restate it")
  refuse(function(x) { x$equations$TF24.assimilation$reads <- list("leafy:carbon"); x },
         "nor a link to a package under links")
  refuse(function(x) { x$links$phylloptim <- ""; x }, "needs a description")
  refuse(function(x) { x$links$phylloptim$models <- list("NOPE"); x },
         "needs the models it serves")
  # A model the package is not declared for may not link to it.
  refuse(function(x) { x$parameters$FF16.lma$link <- "phylloptim:vcmax_25"; x },
         "linked only for TF24")
  # A linked parameter states its own provenance.
  refuse(function(x) { x$parameters$TF24.stem_P50$source <- NULL; x }, "lacks source")
  refuse(function(x) { x$equations$TF24.light_chain$writes <- list("leafy:PPFD"); x },
         "a link to no package under links")
})

test_that("every link written in free text resolves too", {
  skip_if_not_installed("phylloptim")
  # Notes and reasons name phylloptim's ids in prose; a rename there would
  # otherwise leave the prose pointing at nothing.
  reg <- model_registry()
  target <- yaml::read_yaml(system.file("model-registry.yml", package = "phylloptim"))
  text <- c(unlist(reg$raw$parameters), unlist(reg$raw$equations),
            unlist(reg$raw$not_parameters))
  ids <- unique(unlist(regmatches(text, gregexpr("(?<![:A-Za-z])phylloptim:[A-Za-z_][A-Za-z0-9_]*",
                                                 text, perl = TRUE))))
  ids <- sub("^phylloptim:", "", ids)
  expect_gt(length(ids), 0L)
  missing <- setdiff(ids, c(names(target$parameters), names(target$equations)))
  expect_true(length(missing) == 0L,
              info = paste("not in phylloptim's registry:", toString(missing)))
})

test_that("a phylloptim other than the pinned one is warned about", {
  skip_if_not_installed("phylloptim")
  local_mocked_bindings(registry_pinned_version = function(package) "0.0.1")
  expect_warning(model_registry(), "plant is built against phylloptim 0.0.1")
})
