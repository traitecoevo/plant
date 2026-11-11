#' Check performance on current system using package bench
#'
#' @param strategy_types A list of name strategy types to be tests
#' @param iterations The number of iterations to be run
#' @details For branch-to-branch performance comparisons, run `make` before
#' `devtools::load_all(quiet = TRUE)` so compiled code is rebuilt with the
#' package's intended optimization flags.
#'
#' @return A dataframe of results
#' @export

run_plant_benchmarks <- function(strategy_types = list(FF16 = FF16_Model),
                                 iterations = 1) {
  f_scm <- function(x) {
    p0 <- default_parameters(x, quiet = TRUE)
    p <- add_strategies(p0, trait_matrix(0.0825, "lma"))
    run_scm(p)
    invisible(NULL)
  }

  f_build_schedule <- function(x) {
    p <- default_parameters(x, quiet = TRUE)
    p$strategies <- list(strategy_types[[x]]())
    p$birth_rate <- 0.1
    run_scm(p, refine_schedule = TRUE)
    invisible(NULL)
  }

  message("Running benchmarks via `run_plant_benchmarks`")
  strategy <- names(strategy_types)

  exprs <- list(
    scm = quote(f_scm(strategy)),
    build_schedule = quote(f_build_schedule(strategy))
  )
  
  bench::press(strategy = strategy,
               {
                 do.call(
                   bench::mark,
                   c(
                     list(
                       check = FALSE,
                       # We're not expecting different results to be equivalent
                       iterations = iterations
                     ),
                     exprs
                   )
                 )
               })
}

