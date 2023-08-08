## TODO: Test introduce_new_nodes(vector<double>)

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  context(sprintf("Patch-%s",x))

  # initialise birth rate per species
  s <- strategy_types[[x]]()
  s$birth_rate_y <- 1
  s$is_variable_birth_rate <- FALSE
  
  e <- environment_types[[x]]
  plant <- Individual(x, e)(s)
  node <- Node(x, e)(s)
  
  p <- Parameters(x, e)(strategies=list(s),
                        patch_type = 'meta-population')
  
  env <- make_environment(x)
  
  ctrl <- Control()
  patch <- Patch(x, e)(p, env, ctrl)
  cmp <- Node(x, e)(p$strategies[[1]])

  test_that(sprintf("Basics %s", x), {
    ## TODO: This is something that needs validating: the birth_rate and
    expect_equal(patch$size, 1)
    expect_identical(patch$height_max, cmp$height)
    expect_equal(patch$parameters, p)
    
    expect_is(patch$environment, c(paste0(x, "_Environment"), "R6"))
    expect_identical(patch$environment$time, 0.0)

    expect_equal(length(patch$species), 1)
    expect_is(patch$species[[1]], sprintf("Species<%s,%s>",x,e))
    
    # with no nodes, we only expect env vars
    env_size <- env$ode_size
    #fails here
    env_state <- patch$ode_state
    env_rates <- patch$ode_rates
    expect_equal(patch$ode_size, env_size)
    
    # either 0 or numeric(0)
    if(x %in% c("FF16", "FF16r", "K93")) {
      expect_equal(patch$ode_state, numeric(0))
      expect_equal(patch$ode_rates, numeric(0))
    }
    expect_identical(patch$ode_state, env_state)
    expect_identical(patch$ode_rates, env_rates)
    
    ## Empty environment:
    patch$compute_environment()
    expect_identical(patch$compute_competition(0), 0)
    
    expect_error(patch$introduce_new_node(0), "Invalid value")
    expect_error(patch$introduce_new_node(2), "out of bounds")
    
    # introduce a node and expect different results
    node_size <- Node(x, e)(s)$ode_size
    ode_size = node_size + env_size
    patch$introduce_new_node(1)
    expect_equal(patch$node_ode_size, node_size)
    expect_equal(patch$ode_size, ode_size)
    if (x == "FF16") {
      expect_equal(patch$node_ode_size, 7)
      expect_equal(patch$ode_size, 7)
    }
    ## Then pull this out:
    cmp$compute_initial_conditions(patch$environment, patch$pr_survival(0.0), 
                                   patch$species[[1]]$extrinsic_drivers$evaluate("birth_rate", 0))
     
    ode_state <- c(cmp$ode_state, env_state)
    ode_rates <- c(cmp$ode_rates, env_rates)
    expect_identical(patch$ode_state, ode_state)
    expect_identical(patch$ode_rates, ode_rates)
    if (x == "FF16") {
      expect_equal(ode_state, c(0.3441947, 0.009159, 0, 0, 0, 0, 1.08695), tolerance = 1e-4)
      expect_equal(ode_rates, c(0.3341652, 0.01000000, 0, 5.1781e-09, 9.60270e-07, 0, -0.78726), tolerance = 1e-4)
    }
    y <- patch$ode_state
    patch$set_ode_state(y, 0)
    expect_identical(patch$ode_state, y)
    
    ## NOTE: These should be identical, but are merely equal...
    expect_equal(patch$derivs(y, 0), ode_rates)
    
    ## solver <- solver_from_ode_target(patch, p$control$ode_control)
    ## solver$step()
    ## patch$introduce_new_node(1)
    ## expect_equal(patch$ode_size,
    ##             cmp$ode_size * patch$n_individuals)
    
    patch$reset()
    expect_equal(patch$ode_size, env_size)
    expect_identical(patch$environment$time, 0.0)
    
    t <- patch$environment$time # do via environment only?
    
    ## patch$introduce_new_node(1)
    ## h <- patch$state("height")[[1]]
    ## while (patch$time < 25) {
    ##   solver$step()
    ##   t <- c(t, patch$time)
    ##   h <- c(h, patch$state("height")[[1]])
    ## }
    
    ## TODO: This is not really a test, but we need to look at this and
    ## see if it makes any sense at all.
    ## if (interactive()) {
    ##   plot(t, h, type="l")
    ##   plot(patch$environment$environment_interpolator$xy, type="l")
    ## }
    
    ## patch$reset()
    ## patch$introduce_new_node(1)
    ## solver <- solver_from_ode_target(patch, p$control$ode_control)
    
    ## tt <- seq(0, 25, length.out=26)
    ## hh <- patch$state("height")[[1]]
    ## for (ti in tt[-1]) {
    ##   solver$advance(ti)
    ##   hh <- c(hh, patch$state("height")[[1]])
    ## }
    
    ## if (interactive()) {
    ##   plot(t, h, type="l")
    ##   points(tt, hh)
    ## }
    
    ## expect_equal(hh, spline(t, h, xout=tt)$y, tolerance=1e-7)
    
    ## test_that("OK at end of sequence", {
    ##   expect_identical(patch$time, tt[[length(tt)]])
    ##   solver$advance(tt[[length(tt)]])
    ##   expect_identical(patch$time, tt[[length(tt)]])
    ##   expect_error(solver$advance(tt[[length(tt)]] - 1e-8))
    ## })
    
    ## test_that("State get/set works", {
    ##   patch$reset()
    ##   patch$introduce_new_node(1)
    ##   ode.control <- p$control$ode_control
    ##   ode.control$set_parameters(list(step_size_min = 1e-4))
    ##   solver <- solver_from_ode_target(patch, ode.control)
    ##   while (patch$time < 5) {
    ##     solver$step()
    ##     if (patch$time > patch$n_individuals) {
    ##       patch$introduce_new_node(1)
    ##       solver <- solver_from_ode_target(patch, ode.control)
    ##     }
    ##   }
    ##   patch$compute_rates() # require because we just added offspring
    ##   state <- patch$state
    
    ##   patch2 <- new(PatchNodeTop, patch$parameters)
    ##   expect_error(patch2$state <- state)
    ##   patch2$force_state(state)
    ##   expect_identical(patch2$state, state)
    
    ##   ## Check some things that depend on state make sense:
    ##   expect_identical(patch2$environment$environment_interpolator$xy)
    ##   expect_identical(patch2$time, patch$time)
    ##   expect_identical(patch2$ode_state, patch$ode_state)
    ##   expect_identical(patch2$height, patch$state("height"))
    ##   expect_identical(patch2$ode_rates, patch$ode_rates)
    ## })
  })
  
  test_that("Weibull Disturbance as default", {
    expect_equal(patch$time, 0.0)
    expect_equal(patch$pr_survival(patch$time), 1.0)
    
    expect_equal(patch$disturbance_mean_interval(), 30)
    disturbance <- Weibull_Disturbance_Regime(105.32)
    
    patch$set_time(10)
    expect_equal(patch$time, 10)
    
    expect_identical(patch$pr_survival(10), 
                     disturbance$pr_survival(10))

    # This is how we'd like it but Rcpp wouldn't handle a disturbance pointer
    #expect_is(patch$disturbance_regime, "Disturbance")
  })
  
  test_that("No Disturbance for fixed-time patches", {
    p$patch_type <- "fixed"
    env <- make_environment(x)
    ctrl <- scm_base_control()

    patch <- Patch(x, e)(p, env, ctrl)
    
    expect_identical(patch$time, 0.0)
    expect_identical(patch$pr_survival(patch$time), 1.0)
    
    expect_true(is.na(patch$disturbance_mean_interval()))
    disturbance <- No_Disturbance()
    
    patch$set_time(10)
    expect_equal(patch$time, 10)
    
    expect_identical(patch$pr_survival(10), 1)
    
    expect_identical(patch$pr_survival(10), 
                     disturbance$pr_survival(10))
    
    # This is how we'd like it but Rcpp wouldn't handle a disturbance pointer
    #expect_is(patch$disturbance_regime, "Disturbance")
  })
}
