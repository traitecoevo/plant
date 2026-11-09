## Deprecated entry points to the trade-off functions, kept for one minor
## release. Each returns plant's frozen testing trade-offs (the ones
## default_parameters() attaches) with a deprecation warning. Remove them, and
## this file, at the next minor release.

deprecate_hyperpar <- function(old, new) {
  .Deprecated(msg = sprintf(paste0(
    "%s is deprecated and will be removed in the next minor release of plant. ",
    "Trade-offs now travel with a Parameters object (see ?default_parameters): ",
    "use %s. The recommended trade-off functions and curated species and site ",
    "parameter sets are in the fieldtrial package."), old, new))
}

##' Deprecated: the trade-off (hyperparameter) functions by name
##'
##' Trade-offs now travel with a \code{Parameters} object: see
##' \code{\link{default_parameters}}. These names are kept for one minor
##' release and return plant's frozen \emph{testing} trade-offs with a warning.
##' Migration:
##' \itemize{
##'   \item \code{make_FF16_hyperpar(...)} -> \code{default_hyperpar("FF16", ...)}
##'   \item \code{make_TF24_hyperpar(...)}, \code{make_TF24f_hyperpar(...)},
##'     \code{make_TF24t_hyperpar(...)}, \code{make_TF24w_hyperpar(...)} ->
##'     \code{default_hyperpar("TF24", ...)}
##'   \item \code{make_K93_hyperpar(...)} -> \code{NULL} (K93 has no trade-offs)
##'   \item \code{FF16_hyperpar}, \code{TF24_hyperpar}, ... as the
##'     \code{hyperpar =} argument -> drop the argument and build \code{p} with
##'     \code{default_parameters(model)}, or pass \code{param_hyperpar(p)}
##'   \item \code{make_hyperpar(type)} -> \code{function(...) default_hyperpar(type, ...)}
##'   \item \code{hyperpar(type)} -> \code{default_hyperpar(type)}, or
##'     \code{param_hyperpar(p)} for the trade-offs \code{p} carries
##' }
##' @param ... Coefficient overrides, as for \code{\link{default_hyperpar}};
##'   ignored by \code{make_K93_hyperpar}, whose arguments never reached the
##'   model.
##' @param m,s,filter The trade-off contract: a trait matrix, a strategy, and
##'   whether to drop derived columns sitting at the strategy's defaults.
##' @param type A model name, e.g. \code{"FF16"}.
##' @name hyperpar-deprecated
##' @rdname hyperpar-deprecated
##' @keywords internal
NULL

##' @export
##' @rdname hyperpar-deprecated
make_FF16_hyperpar <- function(...) {
  deprecate_hyperpar("make_FF16_hyperpar()", "default_hyperpar(\"FF16\", ...)")
  default_hyperpar("FF16", ...)
}

##' @export
##' @rdname hyperpar-deprecated
make_TF24_hyperpar <- function(...) {
  deprecate_hyperpar("make_TF24_hyperpar()", "default_hyperpar(\"TF24\", ...)")
  default_hyperpar("TF24", ...)
}

##' @export
##' @rdname hyperpar-deprecated
make_TF24f_hyperpar <- function(...) {
  deprecate_hyperpar("make_TF24f_hyperpar()", "default_hyperpar(\"TF24\", ...)")
  default_hyperpar("TF24", ...)
}

##' @export
##' @rdname hyperpar-deprecated
make_TF24t_hyperpar <- function(...) {
  deprecate_hyperpar("make_TF24t_hyperpar()", "default_hyperpar(\"TF24\", ...)")
  default_hyperpar("TF24", ...)
}

##' @export
##' @rdname hyperpar-deprecated
make_TF24w_hyperpar <- function(...) {
  deprecate_hyperpar("make_TF24w_hyperpar()", "default_hyperpar(\"TF24\", ...)")
  default_hyperpar("TF24", ...)
}

##' @export
##' @rdname hyperpar-deprecated
make_K93_hyperpar <- function(...) {
  deprecate_hyperpar("make_K93_hyperpar()",
                     "hyperpar = NULL (K93 has no trade-offs)")
  identity_hyperpar
}

##' @export
##' @rdname hyperpar-deprecated
FF16_hyperpar <- function(m, s, filter = TRUE) {
  deprecate_hyperpar("FF16_hyperpar", "default_parameters(\"FF16\")")
  default_hyperpar("FF16")(m, s, filter)
}
attr(FF16_hyperpar, "input_traits") <- "narea"

##' @export
##' @rdname hyperpar-deprecated
TF24_hyperpar <- function(m, s, filter = TRUE) {
  deprecate_hyperpar("TF24_hyperpar", "default_parameters(\"TF24\")")
  default_hyperpar("TF24")(m, s, filter)
}

##' @export
##' @rdname hyperpar-deprecated
TF24f_hyperpar <- function(m, s, filter = TRUE) {
  deprecate_hyperpar("TF24f_hyperpar", "default_parameters(\"TF24f\")")
  default_hyperpar("TF24")(m, s, filter)
}

##' @export
##' @rdname hyperpar-deprecated
TF24t_hyperpar <- function(m, s, filter = TRUE) {
  deprecate_hyperpar("TF24t_hyperpar", "default_parameters(\"TF24t\")")
  default_hyperpar("TF24")(m, s, filter)
}

##' @export
##' @rdname hyperpar-deprecated
TF24w_hyperpar <- function(m, s, filter = TRUE) {
  deprecate_hyperpar("TF24w_hyperpar", "default_parameters(\"TF24w\")")
  default_hyperpar("TF24")(m, s, filter)
}

##' @export
##' @rdname hyperpar-deprecated
K93_hyperpar <- function(m, s, filter = TRUE) {
  deprecate_hyperpar("K93_hyperpar", "hyperpar = NULL (K93 has no trade-offs)")
  m
}

##' @export
##' @rdname hyperpar-deprecated
make_hyperpar <- function(type) {
  deprecate_hyperpar("make_hyperpar(type)",
                     "function(...) default_hyperpar(type, ...)")
  ## Fail now on an unknown type, as the old switch() did.
  default_hyperpar(type)
  function(...) default_hyperpar(type, ...)
}

##' @export
##' @rdname hyperpar-deprecated
hyperpar <- function(type) {
  deprecate_hyperpar("hyperpar(type)",
                     "default_hyperpar(type), or param_hyperpar(p)")
  default_hyperpar(type)
}
