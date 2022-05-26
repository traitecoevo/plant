## TODO: This is all just ported over from tree1 and needs splitting
## into units.

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  s <- strategy_types[[x]]()
  e <- environment_types[[x]]

  context(sprintf("Node-%s",x))
  test_that("setup, growth rates", {

    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    expect_is(node, sprintf("Node<%s,%s>", x, e))
    expect_is(node$individual, sprintf("Individual<%s,%s>", x, e))

    env <- test_environment(x, 2 * plant$state("height"),
                            light_env=function(x) rep(1, length(x)))

    ## The big unknown is the growth rate gradient calculation; that is,
    ## the derivative d(dh/dt)/dh.
    growth_rate_given_height <- function(height, plant, env) {
      plant$set_state("height", height)
      plant$compute_rates(env)
      plant$rate("height")
    }

    grad_forward <- function(f, x, dx, ...) {
      (f(x + dx, ...) - f(x, ...)) / dx
    }

    plant$compute_rates(env)
    p2 <- Individual(x, e)(s)

    ## First, a quick sanity check that our little function behaves as
    ## expected:
    expect_equal(growth_rate_given_height(plant$state("height"), p2, env),  plant$rate("height"))

    ## With height:
    ctrl <- s$control
    method_args <- list(d=ctrl$node_gradient_eps,
                        eps=ctrl$node_gradient_eps)

    ## With a plant, manually compute the growth rate gradient using
    ## Richarson extrapolation:
    dgdh_richardson <- numDeriv::grad(growth_rate_given_height, plant$state("height"),
                                      plant=p2, env=env,
                                      method.args=method_args)
    ## And also using plain forward differencing:
    dgdh_forward <- grad_forward(growth_rate_given_height, plant$state("height"),
                                 method_args$eps, plant=p2, env=env)

    ## These agree, but not that much:
    expect_equal(dgdh_forward, dgdh_richardson, tolerance=1e-6)

    ## Now, do this with the node:
    dgdh <- node$growth_rate_gradient(env)

    expect_equal(dgdh, dgdh_forward)

    ## Again with Richardson extrapolation:
        node <- Node(x, e)(s)
    node2 <- Node(x, e)(strategy_types[[x]](control=Control(node_gradient_richardson=TRUE)))
    expect_true(node2$individual$strategy$control$node_gradient_richardson)

    ## NOTE: Not sure why this is not identical: it's either a bug
    ## somewhere or an issue due to reusing intervals.
    dgdh2 <- node2$growth_rate_gradient(env)
    ## expect_identical(dgdh2, dgdh_richardson)
    expect_false(identical(dgdh2, dgdh))

    ## p <- node2$individual
    ## p$compute_rates(env)
    ## f <- function(x) {
    ##   growth_rate_given_height(x, p, env)
    ## }
    ## dgdh3 <- test_gradient_richardson(f, p$set_state("height", ), ctrl$node_gradient_eps,
    ##                                   ctrl$node_gradient_richardson_depth)
    ## dgdh3 - dgdh2

    ## This is entirely optional, but kind of nice to see.
    if (interactive()) {
      hh <- seq(plant$state("height") * 0.5, plant$state("height") * 1.5, length.out=101)
      gr <- sapply(hh, growth_rate_given_height, p2, env)
      p2$set_state("height", plant$state("height"))
      h_focus <- plant$state("height")
      g_focus <- growth_rate_given_height(plant$state("height"), p2, env)
      plot(gr ~ hh, xlab="Height", ylab="Growth rate")
      points(g_focus ~ h_focus, col="red", pch=19)
      ## Intercept by solving y = m*x + c for c => (c = y - m * x).
      abline(g_focus - dgdh * plant$state("height"), dgdh)
    }

  })

  ## TODO: Not done yet:
  ##   * Check that the initial conditions are actually correct
  ##   * Check that the rates computed are actually correct
  test_that("ODE interface", {
    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    env <- test_environment(x, 2 * plant$state("height"),
                            light_env=function(x) rep(1, length(x)))

    node$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 1)
    plant$compute_rates(env)

    nms <- c(plant$ode_names,
             "offspring_produced_survival_weighted", "log_density")
    expect_equal(node$ode_size, length(nms))
    expect_equal(node$ode_names, nms)

    ## Mortality is different because that's what Nodes track
    for( v in setdiff(plant$ode_names, "mortality")) {
      expect_equal(node$individual$state(v), plant$state(v))
    }

    ## Set up plant too:
    pr_estab <- plant$establishment_probability(env)

    y <- plant$ode_state
    g <- plant$rate("height")

    ## Ode *values*:
    cmp <- c(plant$internals$states,
             0, # offspring_produced_survival_weighted
             log(pr_estab / g) # log density
             )
    cmp[which(plant$ode_names == 'mortality')] <- -log(pr_estab)
    expect_equal(node$ode_state, cmp)

    expect_identical(node$fecundity, 0.0);

    ## Ode *rates*:
    cmp <- c(plant$internals$rates,
             ## This is different to the approach in tree1?
             plant$rate("fecundity") * exp(-plant$state("mortality")),
             -plant$rate("mortality") - node$growth_rate_gradient(env))


    expect_equal(node$ode_rates, cmp)
  })

  test_that("leaf area calculations", {
    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    env <- test_environment(x, 10,
                            light_env=function(x) rep(1, length(x)))

    h <- node$height

    expect_equal(node$log_density, -Inf) # zero
    expect_equal(exp(node$log_density), 0.0) # zero
    expect_equal(node$competition_effect, 0) # zero density
    node$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 1)

    expect_equal(node$ode_state[[node$ode_size]], node$log_density)
    density <- exp(node$log_density)
    expect_equal(node$competition_effect, plant$compute_competition(0.0) * density)
    expect_equal(node$compute_competition(h / 2), plant$compute_competition(h / 2) * density)

    h <- 8.0
    plant$set_state("height", h)
    v <- node$ode_state
    v[[1]] <- h
    node$ode_state <- v
    expect_identical(plant$state("height"), h)
    expect_identical(node$height, h)

    expect_equal(node$competition_effect, plant$compute_competition(0) * density)
    expect_equal(node$compute_competition(h / 2), plant$compute_competition(h / 2) * density)
  })
}
