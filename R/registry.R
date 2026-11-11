## The model registry: plant's equations and parameters as data, read from the
## installed inst/model-registry.yml. The file's header is the schema; this
## reader enforces it, and tests/testthat/test-registry.R checks the entries
## against the code.

##' The model registry
##'
##' Reads plant's model registry: one record per equation of the
##' size-structured layer plant's models share (and, as they are added, of each
##' model), per parameter those equations read, and per slot (a quantity the
##' shared layer reads and every model supplies).
##' Each record names the C++ symbol that evaluates or holds it. The registry is
##' a record of the model, not a second implementation of it; documents render
##' their equations from it (see \code{\link{registry_equation}}).
##'
##' @param path The registry file; the installed one by default.
##' @return A list of class \code{plant_model_registry} with \code{package},
##'   \code{sections} (a named character vector, in document order),
##'   \code{slots}, \code{parameters} and \code{equations} (data frames, one row
##'   per entry, keyed by \code{id}; \code{reads}, \code{writes} and \code{also}
##'   are list columns), \code{not_parameters} (a named character vector of
##'   reasons), \code{raw} (the parsed YAML) and \code{path}.
##' @export
##' @examples
##' reg <- model_registry()
##' reg$sections
##' registry_equation("density_pde_height", reg)
model_registry <- function(path = system.file("model-registry.yml",
                                              package = "plant")) {
  raw <- registry_read(path)
  structure(list(
    package = raw$package,
    sections = unlist(raw$sections),
    slots = registry_frame(raw$slots,
                           c("symbol", "units", "meaning", "code")),
    parameters = registry_frame(raw$parameters,
                                c("section", "symbol", "units", "class",
                                  "meaning", "source", "placeholder", "code",
                                  "cite", "literal", "field"),
                                list_cols = "aliases"),
    equations = registry_frame(raw$equations,
                               c("section", "title", "latex", "code", "kind",
                                 "fills", "cite", "note", "overstorey"),
                               list_cols = c("reads", "writes", "also")),
    not_parameters = unlist(raw$not_parameters),
    raw = raw,
    path = path),
    class = "plant_model_registry")
}

##' @export
print.plant_model_registry <- function(x, ...) {
  cat(sprintf("<%s model registry: %d sections, %d slots, %d parameters, %d equations>\n",
              x$package, length(x$sections), nrow(x$slots),
              nrow(x$parameters), nrow(x$equations)))
  invisible(x)
}

##' One equation's LaTeX, from the model registry
##'
##' The display a document renders for an equation, so that the page's LaTeX is
##' the registry's. Write it into a page as a display, for example in a chunk
##' with \code{results = "asis"}: \code{cat("$$", registry_equation(id), "$$")}.
##'
##' @param id The equation's id.
##' @param registry A registry from \code{\link{model_registry}}.
##' @return The LaTeX, a single string without surrounding \code{$$}.
##' @export
registry_equation <- function(id, registry = model_registry()) {
  i <- match(id, registry$equations$id)
  if (is.na(i)) {
    stop(sprintf("No equation '%s' in the registry", id), call. = FALSE)
  }
  trimws(registry$equations$latex[[i]])
}

##' Parameter values, from the live model
##'
##' The registry stores no values; this reads each parameter's current default
##' off the model, so a table rendered from the registry shows the number the
##' model runs with. Parameters held on \code{Parameters} are read from
##' \code{scm_base_parameters(model)}; a \code{literal} gives its own value.
##' Parameters held elsewhere are \code{NA}.
##'
##' @param registry A registry from \code{\link{model_registry}}.
##' @param model The model whose defaults to read, e.g. \code{"FF16"}.
##' @return \code{registry$parameters} with a \code{value} column.
##' @export
registry_values <- function(registry = model_registry(), model = "FF16") {
  pars <- registry$parameters
  p <- scm_base_parameters(model)
  field <- ifelse(is.na(pars$field), sub("^.*::", "", pars$code), pars$field)
  on_parameters <- grepl("^plant::Parameters::", pars$code)
  value <- rep(NA_real_, nrow(pars))
  for (i in which(on_parameters)) {
    v <- p[[field[[i]]]]
    if (is.numeric(v) && length(v) == 1L) value[[i]] <- v
  }
  lit <- !is.na(pars$literal)
  value[lit] <- as.numeric(pars$literal[lit])
  pars$value <- value
  pars
}

## Read and validate the registry file against the schema in its header.
registry_read <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    stop("Model registry not found: ", path, call. = FALSE)
  }
  raw <- yaml::read_yaml(path)
  fail <- function(...) stop("Model registry: ", ..., call. = FALSE)
  top <- c("package", "sections", "slots", "parameters", "equations",
           "not_parameters")
  missing <- setdiff(top, names(raw))
  if (length(missing) > 0L) fail("missing top-level keys: ", toString(missing))
  unknown <- setdiff(names(raw), top)
  if (length(unknown) > 0L) fail("unknown top-level keys: ", toString(unknown))

  ids <- c(names(raw$slots), names(raw$parameters), names(raw$equations))
  dup <- unique(ids[duplicated(ids)])
  if (length(dup) > 0L) {
    fail("ids used twice (across slots, parameters and equations): ",
         toString(dup))
  }
  sections <- names(raw$sections)
  # A required field must be present and, except `code` (where ~ marks an
  # identity with no function), not null.
  need <- function(entries, kind, fields) {
    for (id in names(entries)) {
      e <- entries[[id]]
      null <- names(e)[vapply(e, is.null, logical(1))]
      absent <- union(setdiff(fields, names(e)),
                      intersect(setdiff(fields, "code"), null))
      if (length(absent) > 0L) {
        fail(sprintf("%s '%s' lacks %s", kind, id, toString(absent)))
      }
    }
  }
  need(raw$slots, "slot", c("symbol", "units", "meaning", "code"))
  need(raw$parameters, "parameter",
       c("section", "symbol", "units", "class", "meaning", "source",
         "placeholder", "code"))
  need(raw$equations, "equation",
       c("section", "title", "latex", "code", "reads", "writes", "kind"))

  classes <- c("trait", "hyperparameter", "constant", "driver", "setting")
  kinds <- c("closed_form", "identity", "derivative", "integrator")
  cite_ok <- function(x) {
    is.null(x) ||
      (is.character(x) && length(x) == 1L &&
         grepl("^[A-Za-z][A-Za-z-]*-[0-9]{4}[a-z]?$", x))
  }
  for (id in names(raw$parameters)) {
    e <- raw$parameters[[id]]
    if (!e$section %in% sections) fail(sprintf("parameter '%s': unknown section '%s'", id, e$section))
    if (!e$class %in% classes) fail(sprintf("parameter '%s': unknown class '%s'", id, e$class))
    if (!cite_ok(e$cite)) fail(sprintf("parameter '%s': cite '%s' is not Surname-Year", id, e$cite))
    if (!is.null(e$literal) && !(is.numeric(e$literal) && is.finite(e$literal))) {
      fail(sprintf("parameter '%s': literal is not a finite number (write floats as 1.0e-3)", id))
    }
  }
  readable <- c(names(raw$parameters), names(raw$slots))
  for (id in names(raw$equations)) {
    e <- raw$equations[[id]]
    if (!e$section %in% sections) fail(sprintf("equation '%s': unknown section '%s'", id, e$section))
    if (!e$kind %in% kinds) fail(sprintf("equation '%s': unknown kind '%s'", id, e$kind))
    if (!cite_ok(e$cite)) fail(sprintf("equation '%s': cite '%s' is not Surname-Year", id, e$cite))
    bad <- setdiff(unlist(e$reads), readable)
    if (length(bad) > 0L) {
      fail(sprintf("equation '%s' reads %s, which is neither a parameter nor a slot",
                   id, toString(bad)))
    }
    if (!is.null(e$fills) && !e$fills %in% names(raw$slots)) {
      fail(sprintf("equation '%s' fills unknown slot '%s'", id, e$fills))
    }
  }
  for (id in names(raw$not_parameters)) {
    if (!is.character(raw$not_parameters[[id]]) || !nzchar(raw$not_parameters[[id]])) {
      fail(sprintf("not_parameters '%s' needs a reason", id))
    }
  }
  raw
}

## One row per entry, keyed by id. Scalar fields become character columns (NA
## where absent); `list_cols` stay lists of character vectors.
registry_frame <- function(entries, fields, list_cols = character()) {
  ids <- names(entries)
  if (is.null(ids)) ids <- character()
  out <- data.frame(id = ids, stringsAsFactors = FALSE)
  for (f in fields) {
    out[[f]] <- vapply(entries, function(e) {
      v <- e[[f]]
      if (is.null(v)) NA_character_ else as.character(v)
    }, character(1), USE.NAMES = FALSE)
  }
  for (f in list_cols) {
    out[[f]] <- lapply(entries, function(e) as.character(unlist(e[[f]])))
    names(out[[f]]) <- NULL
  }
  out
}
