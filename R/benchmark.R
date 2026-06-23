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

run_plant_benchmarks <- function(strategy_types = list(FF16 = FF16_Strategy),
                                 iterations = 1) {
  f_scm <- function(x) {
    p0 <- scm_base_parameters(x)
    p <- add_strategies(p0, trait_matrix(0.0825, "lma"))
    run_scm(p)
    invisible(NULL)
  }

  f_build_schedule <- function(x) {
    p <- scm_base_parameters(x)
    p$strategies <- list(strategy_types[[x]]())
    p$birth_rate <- 0.1
    run_scm(p, refine_schedule = TRUE)
    invisible(NULL)
  }

  f_mutant <- function(x) {
    p0 <- scm_base_parameters(x)
    p_resident <- add_strategies(p0, trait_matrix(0.0825, "lma"))

    ctrl <- Control()
    ctrl$save_RK45_cache <- TRUE

    scm <- run_scm(p_resident, ctrl = ctrl)

    # One additional mutant strategy around the resident trait value.
    p_mutant <- add_strategies(p_resident, trait_matrix(0.09, "lma"), birth_rate = 1)
    scm$run_mutant(p_mutant)
    invisible(NULL)
  }

  message("Running benchmarks via `run_plant_benchmarks`")
  strategy <- names(strategy_types)

  exprs <- list(
    scm = quote(f_scm(strategy)),
    build_schedule = quote(f_build_schedule(strategy)),
    mutant <- quote(f_mutant(strategy))
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

# Evaluate overheads of having multiple soil layers in TF24, to ensure that
# our approach to computing individual resource consumption is scalable.
run_resource_consumption_benchmarks <- function(its = 10) {
  
  f_scm <- function(layers) {
    p0 <- scm_base_parameters("TF24")
    p0$max_patch_lifetime = 10
    
    p1 <- add_strategies(p0, trait_matrix(0.0825, "lma"))
    
    env <- Environment("TF24")
    env$set_soil_number_of_depths(layers)

    ctrl <- Control()
    out <- run_scm(p1, env, ctrl)
    invisible(NULL)
  }
  
  message("Running resource consumption benchmarks`")
  soil_layers <- c(1, 10, 50, 100)
  bench::press(soil_layers = soil_layers,
               {
                 bench::mark(
                   check = FALSE,
                   # We're not expecting different results to be equivalent
                   iterations = its,
                   scm = f_scm(soil_layers),
                 )
               })
}
