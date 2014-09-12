##' @export
assembler_sample_positive <- function(community0, n_sample=1L,
                                      seed_rain_eps=1e-3,
                                      approximate_type="naive",
                                      ...) {
  births_sys <- make_births_sample_positive(n_sample, approximate_type)
  deaths_sys <- make_deaths_stochastic_naive(seed_rain_eps)
  assembler(community0, births_sys, deaths_sys, ...)
}

make_births_sample_positive <- function(n, approximate_type="naive") {
  approximate_type <- match.arg(approximate_type, c("naive", "gp"))
  function(sys) {
    if (is.null(sys$bounds)) {
      ## No viable region
      ret <- matrix(nrow=0, ncol=1)
      colnames(ret) <- sys$trait_names
    } else {
      sys$add_approximate_landscape(approximate_type)
      f <- sys$landscape_approximate$predict
      ret <- cbind(rejection_sample(n, f, sys$bounds))
      colnames(ret) <- sys$trait_names
      attr(ret, "done") <- nrow(ret) == 0
    }
    ret
  }
}
