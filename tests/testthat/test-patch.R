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
  env <- Environment(x)
  
  ctrl <- Control()
  patch <- Patch(x, e)(p, env, ctrl)
  cmp <- Node(x, e)(p$strategies[[1]])

  test_that(sprintf("Basics %s", x), {
    ## Default birth rate (see #333): a fresh strategy should be a constant
    ## rate of 1 with the variable-rate (spline) path switched off.
    s_default <- strategy_types[[x]]()
    expect_equal(s_default$birth_rate_y, 1)
    expect_false(s_default$is_variable_birth_rate)

    expect_equal(patch$size, 1)
    expect_identical(patch$height_max, cmp$height)
    expect_equal(patch$parameters, p)
    expect_equal(patch$get_area, 1.0)
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
    if(x %in% c("FF16", "K93")) {
      expect_equal(patch$ode_state, numeric(0))
      expect_equal(patch$ode_rates, numeric(0))
    }
    if(x %in% c("TF24")) {
      length_odes <- env$get_soil_number_of_depths()
      soil_moist_inits <- c(rep(env$soil_moist_sat, length_odes)/2, rep(0,4))
      expect_equal(patch$ode_state, soil_moist_inits)
      expect_equal(patch$ode_rates, c(3.312786717, rep(0,4), 1.000000000, 0.996093750, 0.002257735, 0.000000000), tol = 1e-4)
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

    patch$reset()
    expect_equal(patch$ode_size, env_size)
    expect_identical(patch$environment$time, 0.0)
  })

  test_that("change patch size", {
  
    ctrl <- Control()
    e <- environment_types[[x]]
    env <- Environment(x)

    # TF24 is much costlier per step, so use a shorter horizon -- but with a low
    # height-at-maturity (hmat = 5) so plants actually reproduce within it. With
    # the default hmat over so short a patch the TF24 reproductive ratios
    # underflow to ~1e-15 and the checks below would merely compare zeros;
    # hmat = 5 lifts them to O(10) so this is a genuine demographic check.
    # FF16/K93 reproduce fine at their default hmat over the longer horizon.
    if (x == "TF24") {
      max_patch_lifetime <- 4
      mk_strategy <- function()
        add_strategies(scm_base_parameters("TF24"),
                       trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                       hyperpar = TF24_hyperpar, birth_rate = 1)$strategies[[1]]
    } else {
      max_patch_lifetime <- 30
      mk_strategy <- function() strategy_types[[x]]()
    }

    mk_pars <- function(area) Parameters(x, e)(strategies = list(mk_strategy()),
                                               patch_area = area,
                                               max_patch_lifetime = max_patch_lifetime)

    p2  <- mk_pars(2)
    p10 <- mk_pars(10)
    expect_equal(p2$patch_area, 2)
    expect_equal(Patch(x, e)(p2, env, ctrl)$get_area, 2)
    expect_equal(p10$patch_area, 10)
    expect_equal(Patch(x, e)(p10, env, ctrl)$get_area, 10)

    # Same birth rate: the larger patch has less competition, so a higher
    # net reproductive ratio.
    p2$strategies[[1]]$birth_rate_y <- 1
    p10$strategies[[1]]$birth_rate_y <- 1
    expect_warning(scm2   <- run_scm(p2))
    expect_warning(scm10a <- run_scm(p10))
    expect_gt(scm10a$net_reproduction_ratios, scm2$net_reproduction_ratios)

    # Birth rate scaled x5 with the x5 patch area gives the same seed input per
    # unit area, hence the same reproductive ratio as the small patch. (p2 is
    # unchanged from above, so it is not re-run.)
    p10b <- p10
    p10b$strategies[[1]]$birth_rate_y <- 5
    expect_warning(scm10b <- run_scm(p10b))
    expect_equal(scm2$net_reproduction_ratios,
                 scm10b$net_reproduction_ratios, tol = 1e-4)
    expect_gt(scm10a$net_reproduction_ratios, scm10b$net_reproduction_ratios)
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
    env <- Environment(x)
    ctrl <- Control()

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
