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
    ## numDeriv halves its step from d * height four times, so for TF24 it
    ## starts at d = 1e-2 to stay above the leaf solve's noise floor (below).
    rich_args <- if (grepl("TF24", x)) list(d = 1e-2, eps = 1e-2) else method_args
    dgdh_richardson <- numDeriv::grad(growth_rate_given_height, plant$state("height"),
                                      plant=p2, env=env,
                                      method.args=rich_args)
    ## And also using plain forward and backward differencing:
    dgdh_forward <- grad_forward(growth_rate_given_height, plant$state("height"),
                                 method_args$eps, plant=p2, env=env)
    dgdh_backward <- grad_backward(growth_rate_given_height, plant$state("height"),
                                   method_args$eps, plant=p2, env=env)

    ## These agree, but not that much. TF24's growth rate is a smooth-positive
    ## part of net production times the reserve gate (#517), so it carries more
    ## curvature than FF16/K93 and the O(h) forward difference departs from the
    ## Richardson estimate a little more -- allow a looser tolerance there.
    ## TF24's leaf solve also leaves ~1e-11 of noise in the growth rate at the
    ## default stem curve (#654), and at eps = 1e-6 that alone moves the
    ## forward difference by ~4e-5 relative.
    tol_grad <- if (grepl("TF24", x)) 1e-4 else 1e-6
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

  test_that("ODE interface", {
    plant <- Individual(x, e)(s)
    ## The boundary node of a fresh patch, seeded at a birth rate of one and a
    ## patch survival of one, in the patch's environment.
    patch <- one_species_patch(x, s)
    node <- patch$species[[1]]$new_node
    env <- patch$environment

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
    ## The boundary condition differs by coordinate -- a density in height is
    ## divided by the growth rate, a density in birth date is not -- so read the
    ## coordinate off the node rather than assuming one. TF24 resolves Control's
    ## "auto" to birth date; FF16 and K93 to height.
    ## On the birth-date coordinate a node is born per seed that establishes, so
    ## it starts at the birth rate (one here) with no hazard paid; its species
    ## supplies the establishment over its cell and folds it into the state when
    ## the next node arrives.
    birth_date <- node$density_in_birth_date
    cmp <- c(plant$internals$states,
             0, # offspring_produced_survival_weighted
             if (birth_date) 0 else log(pr_estab / g) # log density
             )
    cmp[which(plant$ode_names == 'mortality')] <-
      if (birth_date) 0 else -log(pr_estab)
    expect_equal(node$ode_state, cmp)

    expect_identical(node$fecundity, 0.0);

    ## Ode *rates*:
    ## Only a density in height compresses as the spacing between neighbouring
    ## sizes changes, so only that coordinate carries the gradient term.
    cmp <- c(plant$internals$rates,
             ## This is different to the approach in tree1?
             plant$rate("fecundity") * exp(-plant$state("mortality")),
             -plant$rate("mortality") -
               if (birth_date) 0 else node$growth_rate_gradient(env))


    expect_equal(node$ode_rates, cmp)
  })

  test_that("leaf area calculations", {
    plant <- Individual(x, e)(s)

    ## A node no patch has seeded carries no density.
    unseeded <- Node(x, e)(s)
    expect_equal(unseeded$log_density, -Inf) # zero
    expect_equal(exp(unseeded$log_density), 0.0) # zero
    expect_equal(unseeded$compute_competition(0.0), 0) # zero density

    patch <- one_species_patch(x, s)
    node <- patch$species[[1]]$new_node
    h <- node$height

    ## The ODE state carries the density in its last entry. On the birth-date
    ## coordinate the density the node contributes also counts the establishment
    ## its species supplies; on the height coordinate the two are the same.
    carried <- node$ode_state[[node$ode_size]]
    supplied <- if (node$density_in_birth_date)
      log(node$individual$establishment_probability(patch$environment)) else 0
    expect_equal(node$log_density, carried + supplied)
    density <- exp(node$log_density)
    expect_equal(node$compute_competition(0.0), plant$compute_competition(0.0) * density)
    expect_equal(node$compute_competition(h / 2), plant$compute_competition(h / 2) * density)

    ## Introduced and grown to 8 m through the patch's state: the height reaches
    ## the node's individual, and the density it was born with is kept.
    h <- 8.0
    plant$set_state("height", h)
    patch$introduce_new_node(1, 0)
    set_node_states(patch, "height", h)
    node <- patch$species[[1]]$nodes[[1]]
    expect_identical(node$height, h)
    expect_identical(node$individual$state("height"), h)
    expect_equal(exp(node$log_density), density)

    expect_equal(node$compute_competition(0.0), plant$compute_competition(0) * density)
    expect_equal(node$compute_competition(h / 2), plant$compute_competition(h / 2) * density)
  })
}

for (x in names(strategy_types)) {
  e <- environment_types[[x]]

  test_that(sprintf("density in birth date (%s)", x), {
    s <- strategy_types[[x]]()
    s$birth_rate_y <- 2

    patch <- one_species_patch(x, s, coordinate = "birth_date")
    node <- patch$species[[1]]$new_node
    expect_true(node$density_in_birth_date)

    ## No division by the growth rate at the boundary, and no establishment in
    ## the state: a node is born per seed that establishes, at the birth rate,
    ## with no hazard paid. Its species supplies the establishment over its
    ## cell, which the density the node contributes counts.
    carried <- node$ode_state[[node$ode_size]]
    expect_equal(carried, log(2))
    expect_equal(node$log_density - carried,
                 log(node$individual$establishment_probability(patch$environment)))
    expect_identical(node$individual$state("mortality"), 0)

    ## Nothing moves an individual along the birth-date axis, so the only
    ## term left is mortality.
    rates <- node$ode_rates
    expect_equal(rates[[length(rates)]],
                 -node$individual$rate("mortality"))

    ## And the height coordinate keeps its compression term. Named explicitly:
    ## TF24 resolves "auto" to birth date, so a default Control here would test
    ## the same coordinate as the block above.
    patch2 <- one_species_patch(x, s, coordinate = "height")
    node2 <- patch2$species[[1]]$new_node
    expect_false(node2$density_in_birth_date)
    rates2 <- node2$ode_rates
    expect_equal(rates2[[length(rates2)]],
                 -node2$individual$rate("mortality") -
                   node2$growth_rate_gradient(patch2$environment))
  })
}

## A node's seed output accrues at its fecundity rate times its own survival and
## the patch's survival since its birth: the patch survival now over the patch
## survival it was born at. Checked over a run, whose nodes span a range of
## survivals at birth, so a node that read another's -- or none -- would show.
test_that("the seed-output rate divides by the patch survival at birth", {
  p0 <- testing_parameters("FF16")
  p0$max_patch_lifetime <- 20
  p <- add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
  solver <- DeterministicPatchRunner("FF16", "FF16_Env")(
    p, Environment("FF16"), empty_events(), Control())
  solver$run()

  patch <- solver$patch
  sp <- patch$species[[1]]
  at_birth <- sp$pr_patch_survival_at_birth
  expect_lt(min(at_birth), 0.01)
  expect_equal(max(at_birth), 1)

  pr_now <- patch$pr_survival(patch$time)
  nodes <- sp$nodes
  i <- match("offspring_produced_survival_weighted", nodes[[1]]$ode_names)
  for (k in seq_along(nodes)) {
    ind <- nodes[[k]]$individual
    expect_equal(nodes[[k]]$ode_rates[[i]],
                 ind$rate("fecundity") * exp(-ind$state("mortality")) *
                   pr_now / at_birth[[k]],
                 tolerance = 1e-12)
  }
})
