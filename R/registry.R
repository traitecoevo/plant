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
##' An entry may link to another package's registry instead of restating it:
##' TF24's leaf parameters are phylloptim's, and its equations read phylloptim's
##' equations. A link is written \code{package:id}. A linked parameter takes its
##' symbol, units, class and meaning from its target when the registry is read,
##' so nothing is copied, and states its own source and placeholder, since its
##' value is plant's. In an equation, a link in \code{reads} is something of the
##' other package's that the equation's code reads or calls, and a link in
##' \code{writes} is an input of the other package's that it seats. Links resolve
##' against the installed package, which is the version plant compiles against
##' (\code{LinkingTo} pins it exactly); a target that is missing there is an
##' error, and an installed version other than the pinned one is a warning,
##' because the other registry's ids are not versioned. Where the package is
##' not installed, links stay unresolved and their borrowed fields are
##' \code{NA}.
##'
##' @param path The registry file; the installed one by default.
##' @return A list of class \code{plant_model_registry} with \code{package},
##'   \code{sections} (a named character vector, in document order),
##'   \code{slots}, \code{parameters} and \code{equations} (data frames, one row
##'   per entry, keyed by \code{id}; \code{reads}, \code{writes} and \code{also}
##'   are list columns; a linked parameter's \code{link} names its target),
##'   \code{links} (one row per link: the entry it is on, its target's package
##'   and id, the target's kind, and whether it resolved), \code{linked} (one
##'   row per package linked to: the installed and the pinned version),
##'   \code{not_parameters} (a named character vector of reasons), \code{raw}
##'   (the parsed YAML) and \code{path}.
##' @export
##' @examples
##' reg <- model_registry()
##' reg$sections
##' registry_equation("density_pde_height", reg)
model_registry <- function(path = system.file("model-registry.yml",
                                              package = "plant")) {
  raw <- registry_read(path)
  resolved <- registry_resolve_links(raw)
  structure(list(
    package = raw$package,
    sections = unlist(raw$sections),
    models = registry_models(raw$models),
    slots = registry_frame(raw$slots,
                           c("symbol", "units", "meaning", "code",
                             "per_model")),
    parameters = registry_frame(resolved$parameters,
                                c("model", "section", "symbol", "units",
                                  "class", "meaning", "source", "placeholder",
                                  "code", "cite", "literal", "field", "fills",
                                  "link", "note"),
                                list_cols = "aliases"),
    equations = registry_frame(raw$equations,
                               c("model", "section", "title", "latex", "code",
                                 "kind", "fills", "cite", "note",
                                 "overstorey"),
                               list_cols = c("reads", "writes", "also")),
    links = resolved$links,
    linked = resolved$linked,
    not_parameters = unlist(raw$not_parameters),
    raw = raw,
    path = path),
    class = "plant_model_registry")
}

##' @export
print.plant_model_registry <- function(x, ...) {
  cat(sprintf("<%s model registry: %d models, %d sections, %d slots, %d parameters, %d equations, %d links>\n",
              x$package, nrow(x$models), length(x$sections), nrow(x$slots),
              nrow(x$parameters), nrow(x$equations), nrow(x$links)))
  invisible(x)
}

##' One equation's LaTeX, from the model registry
##'
##' The display a document renders for an equation, so that the page's LaTeX is
##' the registry's. Write it into a page as a display, for example in a chunk
##' with \code{results = "asis"}: \code{cat("$$", registry_equation(id), "$$")}.
##'
##' @param id The equation's id, or a link (\code{"phylloptim:carbon"}) to an
##'   equation in a linked package's registry, whose LaTeX is returned.
##' @param registry A registry from \code{\link{model_registry}}.
##' @return The LaTeX, a single string without surrounding \code{$$}.
##' @export
registry_equation <- function(id, registry = model_registry()) {
  link <- registry_parse_link(id)
  if (!is.null(link)) {
    if (!link[["package"]] %in% registry$linked$package) {
      stop(sprintf("The registry links to no package '%s'", link[["package"]]),
           call. = FALSE)
    }
    target <- registry_target(link[["package"]])
    if (is.null(target)) {
      stop(sprintf("Cannot resolve '%s': %s is not installed", id,
                   link[["package"]]), call. = FALSE)
    }
    e <- target$equations[[link[["id"]]]]
    if (is.null(e)) {
      stop(sprintf("No equation '%s' in %s's registry", link[["id"]],
                   link[["package"]]), call. = FALSE)
    }
    return(trimws(e$latex))
  }
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
##' \code{base_parameters(model)}, and a model's own parameters (those of
##' \code{model} and of any model it extends) from its default strategy's
##' \code{pars}, or from the strategy itself for a subclass's traits held
##' beside \code{pars}; a \code{literal} gives its own value. Other models'
##' parameters, those held on the environment, and any not held as a single
##' number, are \code{NA}.
##'
##' @param registry A registry from \code{\link{model_registry}}.
##' @param model The model whose defaults to read, e.g. \code{"FF16"}.
##' @return \code{registry$parameters} with a \code{value} column.
##' @export
registry_values <- function(registry = model_registry(), model = "FF16") {
  pars <- registry$parameters
  p <- base_parameters(model)
  s <- get(paste0(model, "_Model"), mode = "function")()
  field <- ifelse(is.na(pars$field), sub("^.*::", "", pars$code), pars$field)
  on_parameters <- grepl("^plant::Parameters::", pars$code)
  on_model <- !is.na(pars$model) &
    pars$model %in% registry_lineage(registry$models, model)
  value <- rep(NA_real_, nrow(pars))
  scalar <- function(v) if (is.numeric(v) && length(v) == 1L) v else NA_real_
  for (i in which(on_parameters)) value[[i]] <- scalar(p[[field[[i]]]])
  # A subclass's own traits sit beside pars, at the strategy's top level.
  for (i in which(on_model)) {
    v <- s$pars[[field[[i]]]]
    value[[i]] <- scalar(if (is.null(v)) s[[field[[i]]]] else v)
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
  top <- c("package", "sections", "models", "links", "slots", "parameters",
           "equations", "not_parameters")
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
  need(raw$models, "model", c("title", "header"))
  need(raw$slots, "slot", c("symbol", "units", "meaning", "code"))

  # A link names a declared package and an id in its registry, and is made
  # only by an entry of a model the package is declared for (or one extending
  # it). A linked parameter takes its description from its target and states
  # its own provenance: the value is plant's.
  for (pkg in names(raw$links)) {
    l <- raw$links[[pkg]]
    if (!is.list(l) || !is.character(l$description) || length(l$description) != 1L ||
        !nzchar(l$description)) {
      fail(sprintf("links: '%s' needs a description", pkg))
    }
    bad <- setdiff(unlist(l$models), names(raw$models))
    if (length(l$models) == 0L || length(bad) > 0L) {
      fail(sprintf("links: '%s' needs the models it serves, from models (not %s)",
                   pkg, toString(bad)))
    }
  }
  lineage_of <- function(model) {
    out <- character()
    while (!is.null(model) && !model %in% out) {
      out <- c(out, model)
      model <- raw$models[[model]]$extends
    }
    out
  }
  link_ok <- function(x) {
    l <- registry_parse_link(x)
    !is.null(l) && l[["package"]] %in% names(raw$links)
  }
  link_scoped <- function(kind, id, model, x) {
    pkg <- registry_parse_link(x)[["package"]]
    if (length(intersect(lineage_of(model), unlist(raw$links[[pkg]]$models))) == 0L) {
      fail(sprintf("%s '%s' links to %s, which is linked only for %s", kind, id, x,
                   toString(unlist(raw$links[[pkg]]$models))))
    }
  }
  linked <- vapply(raw$parameters, function(e) !is.null(e$link), logical(1))
  for (id in names(raw$parameters)[linked]) {
    e <- raw$parameters[[id]]
    if (!link_ok(e$link)) {
      fail(sprintf("parameter '%s': link '%s' is not package:id for a package under links",
                   id, toString(e$link)))
    }
    link_scoped("parameter", id, e$model, e$link)
    copied <- intersect(registry_borrowed, names(e))
    if (length(copied) > 0L) {
      fail(sprintf("parameter '%s' links to %s, so it takes %s from there and must not restate it",
                   id, e$link, toString(copied)))
    }
  }
  need(raw$parameters[linked], "parameter",
       c("section", "code", "link", "source", "placeholder"))
  need(raw$parameters[!linked], "parameter",
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
  models <- names(raw$models)
  slots <- names(raw$slots)
  for (id in models) {
    m <- raw$models[[id]]
    if (!is.null(m$extends) && !m$extends %in% models) {
      fail(sprintf("model '%s' extends unknown model '%s'", id, m$extends))
    }
    bad <- setdiff(names(m$unfilled), slots)
    if (length(bad) > 0L) {
      fail(sprintf("model '%s' leaves unknown slots unfilled: %s", id, toString(bad)))
    }
    for (u in names(m$unfilled)) {
      if (!is.character(m$unfilled[[u]]) || !nzchar(m$unfilled[[u]])) {
        fail(sprintf("model '%s' leaves '%s' unfilled with no reason", id, u))
      }
    }
    seen <- id
    parent <- m$extends
    while (!is.null(parent)) {
      if (parent %in% seen) fail(sprintf("model '%s' extends itself through %s", id, toString(seen)))
      seen <- c(seen, parent)
      parent <- raw$models[[parent]]$extends
    }
  }
  for (id in slots) {
    pm <- raw$slots[[id]]$per_model
    if (!is.null(pm) && !isTRUE(pm) && !isFALSE(pm)) {
      fail(sprintf("slot '%s': per_model must be true or false", id))
    }
  }
  # A model's entry says so, and its id is `<model>.<name>`; a shared entry
  # has neither. A shared entry may fill only a slot supplied once for every
  # model (per_model: false), and a model's entry only a per-model slot.
  shared_slots <- slots[vapply(raw$slots, function(x) isFALSE(x$per_model), logical(1))]
  owner_ok <- function(kind, id, e) {
    if (!is.null(e$fills) && !e$fills %in% slots) {
      fail(sprintf("%s '%s' fills unknown slot '%s'", kind, id, e$fills))
    }
    if (is.null(e$model)) {
      if (grepl(".", id, fixed = TRUE)) {
        fail(sprintf("%s '%s' is named like a model's entry but has no model", kind, id))
      }
      if (!is.null(e$fills) && !e$fills %in% shared_slots) {
        fail(sprintf("%s '%s' is shared but fills per-model slot '%s'", kind, id, e$fills))
      }
      return(invisible())
    }
    if (!e$model %in% models) fail(sprintf("%s '%s': unknown model '%s'", kind, id, e$model))
    if (!startsWith(id, paste0(e$model, "."))) {
      fail(sprintf("%s '%s' belongs to %s, so its id starts '%s.'", kind, id, e$model, e$model))
    }
    if (!is.null(e$fills) && e$fills %in% shared_slots) {
      fail(sprintf("%s '%s' fills '%s', which every model shares", kind, id, e$fills))
    }
  }
  for (id in names(raw$parameters)) {
    e <- raw$parameters[[id]]
    owner_ok("parameter", id, e)
    if (!e$section %in% sections) fail(sprintf("parameter '%s': unknown section '%s'", id, e$section))
    if (is.null(e$link) && !e$class %in% classes) {
      fail(sprintf("parameter '%s': unknown class '%s'", id, e$class))
    }
    if (!cite_ok(e$cite)) fail(sprintf("parameter '%s': cite '%s' is not Surname-Year", id, e$cite))
    if (!is.null(e$literal) && !(is.numeric(e$literal) && is.finite(e$literal))) {
      fail(sprintf("parameter '%s': literal is not a finite number (write floats as 1.0e-3)", id))
    }
  }
  readable <- c(names(raw$parameters), names(raw$slots))
  for (id in names(raw$equations)) {
    e <- raw$equations[[id]]
    owner_ok("equation", id, e)
    if (!e$section %in% sections) fail(sprintf("equation '%s': unknown section '%s'", id, e$section))
    if (!e$kind %in% kinds) fail(sprintf("equation '%s': unknown kind '%s'", id, e$kind))
    if (!cite_ok(e$cite)) fail(sprintf("equation '%s': cite '%s' is not Surname-Year", id, e$cite))
    reads <- unlist(e$reads)
    bad <- setdiff(reads, readable)
    bad <- bad[!vapply(bad, link_ok, logical(1))]
    if (length(bad) > 0L) {
      fail(sprintf("equation '%s' reads %s, which is neither a parameter nor a slot, nor a link to a package under links",
                   id, toString(bad)))
    }
    # A link in `writes` is an input of the linked package's that this
    # equation's code seats.
    links_written <- Filter(function(w) !is.null(registry_parse_link(w)),
                            unlist(e$writes))
    bad <- links_written[!vapply(links_written, link_ok, logical(1))]
    if (length(bad) > 0L) {
      fail(sprintf("equation '%s' writes %s, a link to no package under links",
                   id, toString(bad)))
    }
    for (x in c(Filter(link_ok, reads), links_written)) {
      link_scoped("equation", id, e$model, x)
    }
    if (!is.null(e$fills) && !e$fills %in% names(raw$slots)) {
      fail(sprintf("equation '%s' fills unknown slot '%s'", id, e$fills))
    }
  }
  for (id in names(raw$not_parameters)) {
    if (grepl(".", id, fixed = TRUE) && !sub("[.].*$", "", id) %in% models) {
      fail(sprintf("not_parameters '%s' names an unknown model", id))
    }
    if (!is.character(raw$not_parameters[[id]]) || !nzchar(raw$not_parameters[[id]])) {
      fail(sprintf("not_parameters '%s' needs a reason", id))
    }
  }
  raw
}

## The fields a linked parameter takes from its target rather than stating.
## Its source, placeholder and citation are plant's: the value comes off
## plant's model, so its provenance is plant's to state.
registry_borrowed <- c("symbol", "units", "class", "meaning")

## "phylloptim:carbon" -> c(package = "phylloptim", id = "carbon"); NULL if x
## is not written as a link.
registry_parse_link <- function(x) {
  if (!is.character(x) || length(x) != 1L) return(NULL)
  m <- regmatches(x, regexec("^([A-Za-z][A-Za-z0-9.]*):([A-Za-z0-9_.]+)$", x))[[1]]
  if (length(m) == 0L) NULL else c(package = m[[2]], id = m[[3]])
}

## A linked package's registry, as parsed YAML, or NULL where the package is not
## installed. Read from its file rather than through its R reader, whose return
## shape is that package's to change; the YAML keys are its documented schema.
registry_target <- function(package) {
  path <- system.file("model-registry.yml", package = package)
  if (!nzchar(path)) return(NULL)
  yaml::read_yaml(path)
}

## The version plant's LinkingTo pins a package to, or NA.
registry_pinned_version <- function(package) {
  linking <- utils::packageDescription("plant", fields = "LinkingTo")
  if (!is.character(linking) || is.na(linking)) return(NA_character_)
  m <- regmatches(linking, regexec(sprintf("\\b%s\\s*\\(==\\s*([0-9.-]+)\\)", package),
                                   linking))[[1]]
  if (length(m) == 0L) NA_character_ else m[[2]]
}

## Resolve every link against the installed target registries: a linked
## parameter's borrowed fields are filled from its target, and every link is
## listed. A target missing from an installed registry is an error.
registry_resolve_links <- function(raw) {
  fail <- function(...) stop("Model registry: ", ..., call. = FALSE)
  packages <- names(raw$links)
  if (is.null(packages)) packages <- character()
  targets <- lapply(packages, registry_target)
  names(targets) <- packages
  installed <- vapply(packages, function(p) {
    if (is.null(targets[[p]])) NA_character_ else as.character(utils::packageVersion(p))
  }, character(1), USE.NAMES = FALSE)
  pinned <- vapply(packages, registry_pinned_version, character(1), USE.NAMES = FALSE)
  moved <- !is.na(installed) & !is.na(pinned) & installed != pinned
  if (any(moved)) {
    warning(sprintf("Model registry: links resolve against %s, but plant is built against %s; the linked ids may have moved",
                    toString(paste(packages[moved], installed[moved])),
                    toString(paste(packages[moved], pinned[moved]))),
            call. = FALSE)
  }
  linked <- data.frame(package = packages,
                       description = vapply(packages, function(p) raw$links[[p]]$description,
                                            character(1), USE.NAMES = FALSE),
                       models = vapply(packages, function(p) toString(unlist(raw$links[[p]]$models)),
                                       character(1), USE.NAMES = FALSE),
                       installed = installed, pinned = pinned,
                       stringsAsFactors = FALSE)

  rows <- list()
  look_up <- function(from, link, want) {
    l <- registry_parse_link(link)
    t <- targets[[l[["package"]]]]
    kind <- NA_character_
    if (!is.null(t)) {
      kind <- if (!is.null(t$parameters[[l[["id"]]]])) "parameter"
              else if (!is.null(t$equations[[l[["id"]]]])) "equation"
              else NA_character_
      if (is.na(kind)) {
        fail(sprintf("'%s' links to %s, which %s's registry does not hold",
                     from, link, l[["package"]]))
      }
      if (!is.null(want) && kind != want) {
        fail(sprintf("'%s' links to %s, which is %s in %s's registry, not a %s",
                     from, link, if (kind == "equation") "an equation" else "a parameter",
                     l[["package"]], want))
      }
    }
    rows[[length(rows) + 1L]] <<- data.frame(
      from = from, package = l[["package"]], target = l[["id"]], kind = kind,
      resolved = !is.null(t), stringsAsFactors = FALSE)
    if (is.null(t)) NULL else t$parameters[[l[["id"]]]]
  }

  parameters <- raw$parameters
  for (id in names(parameters)) {
    link <- parameters[[id]]$link
    if (is.null(link)) next
    target <- look_up(id, link, "parameter")
    for (f in registry_borrowed) parameters[[id]][[f]] <- target[[f]]
  }
  for (id in names(raw$equations)) {
    e <- raw$equations[[id]]
    for (r in c(unlist(e$reads), unlist(e$writes))) {
      if (!is.null(registry_parse_link(r))) look_up(id, r, NULL)
    }
  }
  links <- if (length(rows) == 0L) {
    data.frame(from = character(), package = character(), target = character(),
               kind = character(), resolved = logical(), stringsAsFactors = FALSE)
  } else {
    do.call(rbind, rows)
  }
  list(parameters = parameters, links = links, linked = linked)
}

## The models, one row each: title, header, the model it extends (NA for none)
## and the slots it leaves unfilled, with why.
registry_models <- function(models) {
  out <- registry_frame(models, c("title", "header", "extends"))
  out$unfilled <- lapply(models, function(m) unlist(m$unfilled))
  names(out$unfilled) <- NULL
  out
}

## A model and every model it extends, nearest first.
registry_lineage <- function(models, model) {
  out <- character()
  while (!is.na(model) && !model %in% out) {
    out <- c(out, model)
    i <- match(model, models$id)
    model <- if (is.na(i)) NA_character_ else models$extends[[i]]
  }
  out
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
