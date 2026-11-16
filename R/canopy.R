##' A strategy's crown shape
##'
##' The vertical distribution of a plant's leaf area, as the model computes it
##' through its crown shape (\code{CanopyShape}): \code{Q}, the fraction of the
##' leaf area above relative height \code{u = z / H}, which is the shade the
##' plant casts at \code{u}; and \code{q = -dQ/du}, the leaf area's density in
##' \code{u}. For FF16, TF24 and K93 alike this is the Yokozawa profile set by
##' \code{pars$eta}. K93's profile spreads basal area over stem diameter rather
##' than leaf area over height.
##'
##' Under FF16's \code{"flat-top-soft-box"} shading model \code{Q} is instead
##' the smoothed step of shade the model has the crown cast, and \code{q} its
##' density; \code{eta_c} is still the Yokozawa profile's, which the model's
##' sapwood and the step's centre read, so it is the mean of \code{u} under
##' \code{q} only where the step lies wholly within the crown
##' (\code{eta_c > 0.5}).
##'
##' The shading model is resolved from the strategy's own \code{control}, with
##' each model's default where it is empty. A strategy taken from a run's
##' parameters carries an empty one, so pass the run's control, which
##' \code{tidy()} keeps as \code{$control}.
##'
##' @param strategy A strategy (a \code{FF16_Model}, \code{TF24_Model},
##'   \code{K93_Model}, ...), such as \code{p$strategies[[1]]}.
##' @param u Relative heights, in [0, 1].
##' @param control The run's \code{Control}, or \code{NULL} for the strategy's
##'   own.
##' @return A list: \code{profile}, a tibble of \code{u}, \code{Q} and
##'   \code{q}; and \code{eta_c}, the Yokozawa profile's leaf-area-weighted
##'   mean of \code{u}.
##' @seealso \code{\link{model_size}}
##' @export
##' @examples
##' crown <- canopy_profile(FF16_Model(), u = seq(0, 1, by = 0.25))
##' crown$profile
##' crown$eta_c
canopy_profile <- function(strategy, u = seq(0, 1, length.out = 101),
                           control = NULL) {
  f <- model_function(strategy, "_canopy_profile")
  if (!is.numeric(u) || any(!is.finite(u) | u < 0 | u > 1)) {
    stop("`u` is a relative height, in [0, 1]", call. = FALSE)
  }
  if (!is.null(control)) {
    if (!inherits(control, "Control")) {
      stop("`control` must be a Control", call. = FALSE)
    }
    strategy$control <- control
  }
  out <- f(strategy, as.numeric(u))
  list(profile = tibble::tibble(u = as.numeric(u), Q = out$Q, q = out$q),
       eta_c = out$eta_c)
}

## The per-model C++ export named <model><suffix>, for a strategy's model.
model_function <- function(strategy, suffix) {
  model <- sub("_Model$", "", class(strategy)[[1]])
  f <- get0(paste0(model, suffix), envir = asNamespace("plant"),
            mode = "function", inherits = FALSE)
  if (is.null(f)) {
    stop("Not a plant model's strategy: ", class(strategy)[[1]], call. = FALSE)
  }
  f
}

##' What a model's size variable is
##'
##' Every plant model integrates one size state, named \code{height} in its
##' outputs, and competes through one size. They are not the same quantity in
##' every model: K93's size is stem diameter and its competitive size is basal
##' area. This reads both, with their units, from the model registry, so a
##' caller labels them without knowing the model.
##'
##' @param x A model name (\code{"FF16"}, \code{"TF24f"}, ...), a
##'   \code{Parameters} object, a \code{plant_run}, or a tidied deterministic
##'   run, which carries its parameters as \code{$p}.
##' @param registry A registry from \code{\link{model_registry}}.
##' @return A list of \code{model}, \code{size} and \code{size_units} (the
##'   state named \code{height}), and \code{competitive_size} and
##'   \code{competitive_size_units}. A model that extends another inherits its
##'   entries.
##' @seealso \code{\link{canopy_profile}}, \code{\link{expand_state}}
##' @export
##' @examples
##' model_size("K93")$size
##' model_size(base_parameters("FF16"))$size_units
model_size <- function(x, registry = model_registry()) {
  if (inherits(x, "plant_run")) {
    x <- x$scenario$p
  } else if (is.list(x) && !inherits(x, "Parameters") && !is.null(x$p)) {
    x <- x$p
  }
  model <- if (inherits(x, "Parameters")) {
    extract_RcppR6_template_types(x, "Parameters")[[1]][1]
  } else if (is.character(x) && length(x) == 1L) {
    x
  } else {
    stop("`x` must be a model name, a Parameters object, a plant_run or a ",
         "tidied deterministic run (a tidied stochastic run has no $p)",
         call. = FALSE)
  }
  models <- registry$models
  if (!model %in% models$id) {
    stop("Unknown model: ", model, call. = FALSE)
  }
  fields <- c("size", "size_units", "competitive_size",
              "competitive_size_units")
  out <- list(model = model)
  for (f in fields) {
    for (m in registry_lineage(models, model)) {
      v <- models[[f]][[match(m, models$id)]]
      if (!is.na(v)) {
        out[[f]] <- v
        break
      }
    }
  }
  out
}
