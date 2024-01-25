

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {

  context(sprintf("Individual-%s",x))

  test_that("Reference comparison", {
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    pl <- Individual(x, e)(s)

    expect_is(pl, sprintf("Individual<%s,%s>",x,e))
    # expect_is(pp, sprintf("IndividualPlus<%s>",x))
    # expect_identical(pp$strategy, s)
    expect_identical(pl$strategy, s)

    ## Expected initial conditions
    h0 <- 10
    pl$set_state("height", h0)

    for(v in pl$ode_names[-1]) {
      expect_identical(pl$state(v), 0.0)
      expect_true(is.na(pl$rate(v)))
    }

    expect_identical(pl$state("height"), h0)

    ## Set and get functions behave identically
    pl$set_state("height", h0)

    expect_identical(pl$state("height"), h0)

    m0 <- 5
    pl$set_state("mortality", m0)
    
    expect_identical(pl$state("mortality"), m0)
    
    f0 <- 8
    pl$set_state("fecundity", f0)
    
    expect_identical(pl$state("fecundity"), f0)
    
    expect_error(pl$state("not a trait name"))
    
    ## Compare internals
    vars_pl <- pl$internals
    
    expect_is(vars_pl, "Internals")
    
    # variable_names <- c("competition_effect", "height", "mortality", "fecundity", "area_heartwood", "mass_heartwood")
    # rate_names <- paste0(setdiff(variable_names, "competition_effect"), "_dt")
    
    # expect_true(all(c(variable_names, rate_names) %in% names(vars_pl)))
    
    ## Compute the vital rates and compare them
    env <- test_environment(x, h0)
    light_env <- attr(env, "light_env") # underlying function
    
    pl$compute_rates(env)
    
    # vars_pp <- pp$internals
    vars_pl <- pl$internals
    
    # expect_equal(vars_pl[variable_names], vars_pp[variable_names])
    # expect_equal(vars_pl[rate_names], vars_pp[rate_names])
    
    ## Area_leaf_above
    # for (h in seq(0, h0, length.out=10)) {
    #   expect_identical(pl$compute_competition(h), pp$compute_competition(h))
    # }
    
    ## Germination_probability
    # expect_identical(pl$establishment_probability(env), pp$establishment_probability(env))
    
    ## ode_system
    # expect_identical(pl$ode_size, 5)
    
    # expect_identical(pl$ode_names, ode_names)
    
    ## ode_state
    if(grepl("K93", x))
      expect_equal(pl$ode_state, c(h0, m0, f0))
    else 
      expect_equal(pl$ode_state, c(h0, m0, f0, 0, 0))
    
    expect_equal(pl$strategy_name, x)
  })
  
  test_that("stochastic support", {
    
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    p <- Individual(x, e)(s)
    
    expect_equal(p$state("mortality"), 0.0)
    expect_equal(p$mortality_probability, 0.0)
    
    p$set_state("mortality", pi)
    expect_equal(p$mortality_probability, 1 - exp(-pi))
    
    p$set_state("mortality", 10000000)
    expect_equal(p$mortality_probability, 1)
    
    p$reset_mortality()
    expect_equal(p$state("mortality"), 0.0)
    expect_equal(p$mortality_probability, 0.0)
  })

  test_that("resource_compensation_point", {
  if(x == "FF16") {
    ## R implementation:
    resource_compensation_point_R <- function(x, plant, ...) {
      target <- function(canopy_openness) {
        env <- fixed_environment(x, canopy_openness)
        plant$compute_rates(env)
        plant$aux("net_mass_production_dt")
      }
      
      f1 <- target(1)
      if (f1 < 0.0) {
        NA_real_
      } else {
        uniroot(target, c(0, 1), f.upper=f1, ...)$root
      }
    }
    
    e <- environment_types[[x]]
    p <- Individual(x, e)(strategy_types[[x]]())
    # skip("Comparison no longer evaluate the nesting is too deep")
    expect_equal(p$resource_compensation_point(), resource_compensation_point_R(x, p), tolerance=1e-5)
  } else {
    expect_equal(x, x) # to avoid a warning for empty
  }
  })


  test_that("Maximise individual rate", {
  if(x %in% c("FF16", "FF16r")) {
    
    #set bounds
    bounds = bounds(lma=c(0.01, 3))
    
    opt_res_height_by_size_ode <- optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 3, type = x, size_name = "height", set_state_by_ode = TRUE)

    #check outcome real and positive
    expect_true(opt_res_height_by_size_ode > 0)
    #check trait optima equal to realistic value
    expect_equal(opt_res_height_by_size_ode[1],0.09015998)
    #check height growth rate equal to realistic value
    expect_equal(attr(opt_res_height_by_size_ode, "height_growth_rate"), 1.141283, tolerance = 1e-6)
    
    
    #check that height and size are convergent with same parameters when using set_state_ode = TRUE
    
    opt_res_height_by_height_ode <- optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 3, type = x, size_name = "height", set_state_by_ode = TRUE)
    expect_true(opt_res_height_by_size_ode == opt_res_height_by_height_ode)  
    
    #check that height opt is convergent when comparing set_state_ode = TRUE and set_state_ode = FALSE
    
    opt_res_height_by_size_grow <- optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 3, type = x, size_name = "height", set_state_by_ode = FALSE)
    expect_equal(opt_res_height_by_size_ode, opt_res_height_by_size_grow, 1e-6)
    
    #check optimise_individual_rate_at_size_by_trait for rate other than height
    opt_res_area_heartwood_by_size_grow <- optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 0.1, type = x, size_name = "area_heartwood", set_state_by_ode = FALSE)
    expect_equal(opt_res_area_heartwood_by_size_grow[1], 0.4160424, 1e-6)
    
    #check what happens when unknown rate name is passed in as size name
    expect_error(optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 0.1, type = x, size_name = "mss_heartwood", set_state_by_ode = FALSE), "key not found")
    
    #check that only one size can be passed to function at a time (grow individual to size can take many values as a fector)
    expect_error(optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = c(3,6,9), type = x, size_name = "height", set_state_by_ode = TRUE), "Expecting a single value")
    
    #set bounds - too low for positive growth
    bounds = bounds(lma=c(1e-8, 1e-7))
    
    #optimum trait should be NA and optimum height growth rate 0
    size_res_ode_state_true_no_pos_growth <- optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 3, type = x, size_name = "height", set_state_by_ode = TRUE)
    expect_true(is.na(size_res_ode_state_true_no_pos_growth[1]))
    expect_equal(attr(size_res_ode_state_true_no_pos_growth, "height_growth_rate"), 0)
    
    
    #optimum trait should be NA and optimum height growth rate 0
    size_res_ode_state_false_no_pos_growth <- optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 3, type = x, size_name = "height", set_state_by_ode = TRUE)
    expect_true(is.na(size_res_ode_state_false_no_pos_growth[1]))
    expect_equal(attr(size_res_ode_state_false_no_pos_growth, "height_growth_rate"), 0)
    
    #check when bounds allow positive growth but dont include root (approx 0.09)
    bounds = bounds(lma=c(0.01, 0.07))
    size_res_ode_state_true_no_root = optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 3, type = x, size_name = "height", set_state_by_ode = TRUE)
    expect_equal(0.07, size_res_ode_state_true_no_root[1], 1e-6)
    
    size_res_ode_state_false_no_root = optimise_individual_rate_at_size_by_trait(bounds, log_scale = TRUE, tol = 0.001, size = 3, type = x, size_name = "height", set_state_by_ode = FALSE)
    expect_equal(0.07, size_res_ode_state_true_no_root[1], 1e-6)
    } else {
    expect_equal(x, x) # to avoid a warning for empty
  }
})
}
