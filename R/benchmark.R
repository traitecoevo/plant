#' Check performance on current system using package bench
#'
#' @param strategy_types A list of name strategy types to be tests
#' @param iterations The number of iterations to be run
#'
#' @return A dataframe of results
#' @export

run_plant_benchmarks <- function(strategy_types = list(FF16 = FF16_Strategy,
                                                       # FF16r = FF16r_Strategy,
                                                       FF16w = FF16w_Strategy,
                                                       K93 = K93_Strategy),
                                 iterations = 1) {
  f_scm <- function(x) {
    p0 <- scm_base_parameters(x)
    p <- expand_parameters(trait_matrix(0.0825, "lma"), p0, mutant = FALSE)
    res <- run_scm(p)
  }

  f_build_schedule <- function(x) {
    p <- scm_base_parameters(x)
    p$strategies <- list(strategy_types[[x]]())
    p$birth_rate <- 0.1
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

run_resource_consumption_benchmarks <- function(its = 10) {
  
  f_scm <- function(layers) {
    p0 <- scm_base_parameters("FF16w")
    p0$max_patch_lifetime = 10
    
    p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16w_hyperpar,FALSE)
    
    env <- make_environment("FF16w", 
                            soil_number_of_depths = layers,
                            soil_initial_state = rep(1, layers))
    
    ctrl <- scm_base_control()
    env$set_extrinsic_driver("rainfall", 0:10, 0:10)
    
    out <- run_scm(p1, env, ctrl)
  }
  
  message("Running resource consumption benchmarks`")
  bench::press(soil_layers = c(1, 10, 50, 100),
               {
                 bench::mark(
                   check = FALSE,
                   # We're not expecting different results to be equivalent
                   iterations = its,
                   scm = f_scm(soil_layers),
                 )
               })
}
