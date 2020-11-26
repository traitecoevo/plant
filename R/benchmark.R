#' Check performance on current system using package bench
#'
#' @param strategy_types A list of name strategy types to be tests
#' @param iterations The number of iterations to be run
#'
#' @return A dataframe of results
#' @export

run_plant_benchmarks <- function(strategy_types = list(FF16 = FF16_Strategy,
                                                       FF16r = FF16r_Strategy,
                                                       K93 = K93_Strategy),
                                 iterations = 1) {
  f_scm <- function(x) {
    p0 <- scm_base_parameters(x)
    p <-
      expand_parameters(trait_matrix(0.0825, "lma"), p0, mutant = FALSE)
    if (grepl("K93", x))
      p$k_I <- 1e-3
    res <- run_scm(p)
  }
  
  f_build_schedule <- function(x) {
    p <- scm_base_parameters(x)
    p$strategies <- list(strategy_types[[x]]())
    p$seed_rain <- 0.1
    
    if (grepl("K93", x))
      p$k_I <- 1e-3
    
    p <- build_schedule(p)
  }
  
  message("Running benchmarks via `run_plant_benchmarks`")
  bench::press(strategy = names(strategy_types),
               {
                 bench::mark(
                   check = FALSE,
                   # We're not expecting different results to be equivalent
                   iterations = iterations,
                   scm = f_scm(strategy),
                   build_schedule = f_build_schedule(strategy)
                 )
               })
}
