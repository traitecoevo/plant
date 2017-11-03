## TODO: This is all just ported over from tree1 and needs splitting
## into units.

strategy_types <- get_list_of_strategy_types()

for (x in names(strategy_types)) {

  context(sprintf("Cohort-%s",x))
  test_that("Ported from tree1", {

    s <- strategy_types[[x]]()
    plant <- Plant(x)(s)
    cohort <- Cohort(x)(s)

    expect_is(cohort, sprintf("Cohort<%s>", x))
    expect_is(cohort$plant, sprintf("Plant<%s>", x))

    env <- test_environment(2 * plant$height,
                            light_env=function(x) rep(1, length(x)),
                            seed_rain=1.0)

    ## The big unknown is the growth rate gradient calculation; that is,
    ## the derivative d(dh/dt)/dh.
    growth_rate_given_height <- function(height, plant, env) {
      plant$height <- height
      plant$compute_vars_phys(env)
      plant$internals[["height_dt"]]
    }
    grad_forward <- function(f, x, dx, ...) {
      (f(x + dx, ...) - f(x, ...)) / dx
    }

    plant$compute_vars_phys(env)
    p2 <- PlantPlus(x)(s)

    ## First, a quick sanity check that our little function behaves as
    ## expected:
    expect_equal(growth_rate_given_height(plant$height, p2, env), plant$internals[["height_dt"]])

    ## With height:
    ctrl <- s$control
    method_args <- list(d=ctrl$cohort_gradient_eps,
                        eps=ctrl$cohort_gradient_eps)

    ## With a plant, manually compute the growth rate gradient using
    ## Richarson extrapolation:
    dgdh_richardson <- numDeriv::grad(growth_rate_given_height, plant$height,
                                      plant=p2, env=env,
                                      method.args=method_args)
    ## And also using plain forward differencing:
    dgdh_forward <- grad_forward(growth_rate_given_height, plant$height,
                                 method_args$eps, plant=p2, env=env)

    ## These agree, but not that much:
    expect_equal(dgdh_forward, dgdh_richardson, tolerance=1e-6)

    ## Now, do this with the cohort:
    dgdh <- cohort$growth_rate_gradient(env)

    expect_equal(dgdh, dgdh_forward)

    ## Again with Richardson extrapolation:
        cohort <- Cohort(x)(s)
    cohort2 <- Cohort(x)(strategy_types[[x]](control=Control(cohort_gradient_richardson=TRUE)))
    expect_true(cohort2$plant$strategy$control$cohort_gradient_richardson)

    ## NOTE: Not sure why this is not identical: it's either a bug
    ## somewhere or an issue due to reusing intervals.
    dgdh2 <- cohort2$growth_rate_gradient(env)
    ## expect_identical(dgdh2, dgdh_richardson)
    expect_not_identical(dgdh2, dgdh)

    ## p <- cohort2$plant
    ## p$compute_vars_phys(env)
    ## f <- function(x) {
    ##   growth_rate_given_height(x, p, env)
    ## }
    ## dgdh3 <- test_gradient_richardson(f, p$height, ctrl$cohort_gradient_eps,
    ##                                   ctrl$cohort_gradient_richardson_depth)
    ## dgdh3 - dgdh2

    ## This is entirely optional, but kind of nice to see.
    if (interactive()) {
      hh <- seq(plant$height * 0.5, plant$height * 1.5, length.out=101)
      gr <- sapply(hh, growth_rate_given_height, p2, env)
      p2$height <- plant$height
      h_focus <- plant$height
      g_focus <- growth_rate_given_height(plant$height, p2, env)
      plot(gr ~ hh, xlab="Height", ylab="Growth rate")
      points(g_focus ~ h_focus, col="red", pch=19)
      ## Intercept by solving y = m*x + c for c => (c = y - m * x).
      abline(g_focus - dgdh * plant$height, dgdh)
    }

  })

  ## TODO: Not done yet:
  ##   * Check that the initial conditions are actually correct
  ##   * Check that the rates computed are actually correct
  test_that("ODE interface", {
    s <- strategy_types[[x]]()
    plant <- PlantPlus(x)(s)
    cohort <- Cohort(x)(s)

    env <- test_environment(2 * plant$height,
                            light_env=function(x) rep(1, length(x)),
                            seed_rain=1.0)

    cohort$compute_initial_conditions(env)
    plant$compute_vars_phys(env)

    expect_equal(cohort$ode_size, 6)
    nms <- c("height", "mortality",
             "area_heartwood", "mass_heartwood",
             "seeds_survival_weighted", "log_density")
    expect_equal(cohort$ode_names, nms)

    ## Mortality is different because that's what Cohorts track
    v <- setdiff(names(cohort$plant$internals), "mortality")
    expect_equal(cohort$plant$internals[v], plant$internals[v])

    ## Set up plant too:
    pr_germ <- plant$germination_probability(env)

    y <- plant$ode_state
    g <- plant$internals[["height_dt"]]

    ## Ode *values*:
    cmp <- c(plant$height,
             -log(pr_germ),
             0, # area_heardwood
             0, # mass_heartwood
             0.0, # fecundity
             log(pr_germ * env$seed_rain_dt / g))
    expect_equal(cohort$ode_state, cmp)

    expect_identical(cohort$fecundity, 0.0);

    ## Ode *rates*:
    env$time <- 10
    patch_survival <- env$patch_survival

    rates <- plant$internals

    cmp <- c(rates[["height_dt"]],
             rates[["mortality_dt"]],
             rates[["area_heartwood_dt"]],
             rates[["mass_heartwood_dt"]],
             ## This is different to the approach in tree1?
             rates[["fecundity_dt"]] *
               patch_survival * exp(-cohort$plant$mortality),
             -rates[["mortality_dt"]] - cohort$growth_rate_gradient(env))
    expect_equal(cohort$ode_rates, cmp)
  })

  test_that("leaf area calculations", {
    s <- strategy_types[[x]]()
    plant <- PlantPlus(x)(s)
    cohort <- Cohort(x)(s)

    env <- test_environment(10,
                            light_env=function(x) rep(1, length(x)),
                            seed_rain=1.0)

    h <- cohort$height

    expect_equal(cohort$log_density, -Inf) # zero
    expect_equal(exp(cohort$log_density), 0.0) # zero
    expect_equal(cohort$area_leaf, 0) # zero density
    cohort$compute_initial_conditions(env)
    expect_equal(cohort$ode_state[[6]], cohort$log_density)
    density <- exp(cohort$log_density)
    expect_equal(cohort$area_leaf, plant$area_leaf * density)
    expect_equal(cohort$area_leaf_above(h / 2), plant$area_leaf_above(h / 2) * density)

    h <- 8.0
    plant$height <- h
    v <- cohort$ode_state
    v[[1]] <- h
    cohort$ode_state <- v
    expect_identical(plant$height, h)
    expect_identical(cohort$height, h)

    expect_equal(cohort$area_leaf, plant$area_leaf * density)
    expect_equal(cohort$area_leaf_above(h / 2), plant$area_leaf_above(h / 2) * density)
  })
}
