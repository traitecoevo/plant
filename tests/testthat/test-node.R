strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

for (x in names(strategy_types)) {
  s <- strategy_types[[x]]()
  e <- environment_types[[x]]

  test_that("setup, growth rates", {

    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    expect_inherits(node, sprintf("Node<%s,%s>", x, e))
    expect_inherits(node$individual, sprintf("Individual<%s,%s>", x, e))

    env <- Environment(x)
    env$set_fixed_environment(1.0, 100)

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

    grad_backward <- function(f, x, dx, ...) {
      (f(x, ...) - f(x - dx, ...)) / dx
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
    ## And also using plain forward and backward differencing:
    dgdh_forward <- grad_forward(growth_rate_given_height, plant$state("height"),
                                 method_args$eps, plant=p2, env=env)
    dgdh_backward <- grad_backward(growth_rate_given_height, plant$state("height"),
                                   method_args$eps, plant=p2, env=env)

    ## These agree, but not that much. TF24's growth rate is a smooth-positive
    ## part of net production times the reserve gate (#517), so it carries more
    ## curvature than FF16/K93 and the O(h) forward difference departs from the
    ## Richardson estimate a little more -- allow a looser tolerance there.
    tol_grad <- if (grepl("TF24", x)) 1e-5 else 1e-6
    expect_equal(dgdh_forward, dgdh_richardson, tolerance = tol_grad)

    ## Now, do this with the node. The default control uses backward
    ## differencing (node_gradient_direction = -1), so it matches the
    ## backward difference exactly.
    dgdh <- node$growth_rate_gradient(env)

    expect_equal(dgdh, dgdh_backward)

    ## Again with Richardson extrapolation:
        node <- Node(x, e)(s)
    node2 <- Node(x, e)(strategy_types[[x]](control=Control(node_gradient_richardson=TRUE)))
    expect_true(node2$individual$strategy$control$node_gradient_richardson)

  })

  ## TODO(#482): Not done yet:
  ##   * Check that the initial conditions are actually correct
  ##   * Check that the rates computed are actually correct
  test_that("ODE interface", {
    plant <- Individual(x, e)(s)
    node <- Node(x, e)(s)

    env <- Environment(x)
    env$set_fixed_environment(1.0, 100)

    node$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 1)
    ## Seed the reference plant's strategy-specific initial states the same way
    ## the node does (via compute_initial_conditions), so the two agree for
    ## strategies that carry seeded states -- e.g. TF24's NSC storage pool (#517).
    plant$set_initial_states(env)
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

    env <- Environment(x)
    env$set_fixed_environment(1, 100)

    h <- node$height

    expect_equal(node$log_density, -Inf) # zero
    expect_equal(exp(node$log_density), 0.0) # zero
    expect_equal(node$compute_competition(0.0), 0) # zero density
    node$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 1)

    expect_equal(node$ode_state[[node$ode_size]], node$log_density)
    density <- exp(node$log_density)
    expect_equal(node$compute_competition(0.0), plant$compute_competition(0.0) * density)
    expect_equal(node$compute_competition(h / 2), plant$compute_competition(h / 2) * density)

    h <- 8.0
    plant$set_state("height", h)
    v <- node$ode_state
    v[[1]] <- h
    node$ode_state <- v
    expect_identical(plant$state("height"), h)
    expect_identical(node$height, h)

    expect_equal(node$compute_competition(0.0), plant$compute_competition(0) * density)
    expect_equal(node$compute_competition(h / 2), plant$compute_competition(h / 2) * density)
  })
}

for (x in names(strategy_types)) {
  e <- environment_types[[x]]

  test_that(sprintf("density in birth date (%s)", x), {
    ctrl <- Control()
    ctrl$node_density_in_birth_date <- TRUE
    s <- strategy_types[[x]]()
    s$control <- ctrl

    env <- Environment(x)
    env$set_fixed_environment(1.0, 100)
    env$time <- 3.5

    node <- Node(x, e)(s)
    node$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 2)

    ## No division by the growth rate at the boundary.
    pr_estab <- node$individual$establishment_probability(env)
    expect_equal(node$log_density, log(2 * pr_estab))

    ## Nothing moves an individual along the birth-date axis, so the only
    ## term left is mortality.
    node$compute_rates(env, 1)
    rates <- node$ode_rates
    expect_equal(rates[[length(rates)]],
                 -node$individual$rate("mortality"))

    ## And the height coordinate keeps its compression term.
    s2 <- strategy_types[[x]]()
    node2 <- Node(x, e)(s2)
    node2$compute_initial_conditions(env, pr_patch_survival = 1, birth_rate = 2)
    node2$compute_rates(env, 1)
    rates2 <- node2$ode_rates
    expect_equal(rates2[[length(rates2)]],
                 -node2$individual$rate("mortality") -
                   node2$growth_rate_gradient(env))
  })
}
