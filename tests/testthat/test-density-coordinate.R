# Carrying the SCM's size distribution in birth date rather than in height
# (Control$node_density_in_birth_date).
#
# The two coordinates describe the same population, so for a strategy whose
# growth is a function of size alone they must agree in the limit of a fine
# node schedule. They do not agree at a coarse one, and the gap is almost
# entirely the *height* coordinate's error: see the convergence tests below,
# where the birth-date answer is already converged at the default schedule
# while the height answer is still climbing toward it.
#
# The defect this file exists to catch is a stale or duplicated quadrature
# abscissa: a resource integral taken over the wrong axis, or over an axis with
# zero-width intervals. The first group of tests asserts that directly, on a
# single patch state, naming the defect each one catches and demonstrating (with
# expect_failure) that the assertion really does fail when the defect is
# present. Those cost nothing. The convergence trends that follow are the
# end-to-end check that the invariants compose, through a full run and a shared
# competition profile, into two coordinates that agree; they cost SCM runs, so
# each is carried at the cheapest schedule and patch lifetime that still shows
# its signal.

interleave_schedule <- function(p, n = 1) {
  for (i in seq_len(n)) {
    p$node_schedule_times <- lapply(p$node_schedule_times, function(tt) {
      sort(unique(c(tt, (head(tt, -1) + tail(tt, -1)) / 2)))
    })
  }
  p
}

size_only_parameters <- function(x) {
  p0 <- scm_base_parameters(x)
  if (x == "FF16") {
    add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
  } else {
    add_strategies(p0, trait_matrix(0.059, "b_0"), birth_rate = 1.0)
  }
}

offspring_in_coordinate <- function(p, x, birth_date, step_max = NULL) {
  ctrl <- Control()
  ctrl$node_density_in_birth_date <- birth_date
  if (!is.null(step_max)) {
    ctrl$ode_step_size_max <- step_max
  }
  run_scm(p, Environment(x), ctrl)$offspring_production
}

## A short collected run per (model, coordinate), computed once and shared.
## Every test below only reads these, so memoising them is safe -- and it is
## what makes the invariant tests free: they interrogate a patch state that has
## already been paid for rather than commissioning a run of their own.
short_run <- local({
  cache <- list()
  function(x, birth_date, lifetime = 20) {
    key <- paste(x, birth_date, lifetime, sep = "/")
    if (is.null(cache[[key]])) {
      p0 <- scm_base_parameters(x)
      p0$max_patch_lifetime <- lifetime
      p <- if (x == "FF16") {
        add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
      } else {
        add_strategies(p0, trait_matrix(0.059, "b_0"), birth_rate = 1.0)
      }
      ctrl <- Control()
      ctrl$node_density_in_birth_date <- birth_date
      scm <- SCM(x, environment_type(x))(p, Environment(x), ctrl)
      scm$collect <- TRUE
      scm$run()
      cache[[key]] <<- scm
    }
    cache[[key]]
  }
})

collected_run <- function(birth_date, lifetime = 20) {
  short_run("FF16", birth_date, lifetime)
}

## Sample a handful of steps spanning a run: the two coordinates' grids drift
## apart as the node list ages, so an invariant that holds only at the end is
## not enough.
sampled_steps <- function(scm, n = 4) {
  unique(round(seq(2, length(scm$history), length.out = n)))
}

## ---------------------------------------------------------------------------
## The invariants, asserted directly on one patch state. Each is O(1) in the
## node schedule, and each is paired with the defect that makes it fail.
## ---------------------------------------------------------------------------

## The trapezium rule, written out. Comparing it against
## Species::compute_competition() tests *which grid* the C++ walks; it is only a
## restatement of the C++ if the grid is the same one, which is the point.
trapezium_by_hand <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1))) / 2
}

## The abscissa Species::compute_competition() claims to integrate over:
## introduction times on the birth-date path, negated heights on the height one
## (negated so both increase as the node list is walked from the tallest down).
## The boundary node closes the grid in either case.
quadrature_grid <- function(s, birth_date = s$density_in_birth_date) {
  if (birth_date) {
    c(s$node_times, s$new_node$introduction_time)
  } else {
    -c(s$heights, s$new_node$height)
  }
}

## Each node's *carried* density times its crown area above z: the integrand of
## the competition integral, whichever coordinate is carried.
competition_integrand <- function(s, z) {
  c(vapply(s$nodes, function(nd) nd$compute_competition(z), numeric(1)),
    s$new_node$compute_competition(z))
}

## The two assertions, as helpers, so that the code guarding each invariant is
## literally the code shown to fail when the defect is present.
expect_integrates_over <- function(s, grid, z) {
  expect_equal(s$compute_competition(z),
               trapezium_by_hand(grid, competition_integrand(s, z)),
               tolerance = 1e-10)
}
expect_strictly_increasing <- function(x) {
  expect_true(all(diff(x) > 0))
}

## Defect: the competition integral taken over the wrong axis -- a density per
## unit birth date integrated over height, or the reverse. Caught here
## immediately and by name; the convergence trends at the end of the file catch
## it only slowly and by implication, as a gap that fails to shrink.
test_that("the competition integral runs over the coordinate the density is carried in", {
  for (x in c("FF16", "K93")) {
    for (birth_date in c(FALSE, TRUE)) {
      scm <- short_run(x, birth_date)
      for (k in sampled_steps(scm)) {
        s <- scm$history[[k]]$species[[1]]
        expect_identical(s$density_in_birth_date, birth_date)
        for (z in c(0, s$height_max * 0.3, s$height_max * 0.7)) {
          expect_integrates_over(s, quadrature_grid(s), z)
        }
      }
      ## And the same assertion over the *other* coordinate's grid fails, so it
      ## discriminates between the two axes rather than holding for either.
      ## Measured discrepancy at z = 0: >=20% carrying in birth date, >=0.7%
      ## carrying in height -- seven orders above the tolerance above.
      s <- scm$patch$species[[1]]
      expect_failure(
        expect_integrates_over(s, quadrature_grid(s, !birth_date), 0))
    }
  }
})

## Defect: repeated abscissae. Those span zero width, so the nodes between them
## drop out of the integral -- silently, because a trapezium sum over a grid
## with a zero-width interval is still a valid trapezium sum over that grid. No
## comparison of integral values can see it (the assertion above holds either
## way, as the next test shows), so the grid is checked for strict monotonicity
## in its own right.
test_that("the birth-date quadrature grid is strictly increasing", {
  for (x in c("FF16", "K93")) {
    scm <- short_run(x, TRUE)
    for (k in sampled_steps(scm)) {
      s <- scm$history[[k]]$species[[1]]
      if (s$size < 2) next
      grid <- quadrature_grid(s)
      ## The introduced nodes: no zero-width interval anywhere.
      expect_strictly_increasing(head(grid, -1))
      ## The boundary node closes the grid. Its birth date is the current patch
      ## time, so it is never *behind* the newest node -- but it coincides with
      ## it at the instant of introduction, where a zero-width closing segment
      ## is legitimate and contributes nothing.
      expect_gte(tail(grid, 1), tail(head(grid, -1), 1))
    }
  }
})

## What a repeated birth date actually costs, and why the check above is the one
## that catches it. Introducing the boundary node twice without advancing the
## clock is the only way to produce one, and it is what a patch seeded without
## per-node times does.
test_that("a repeated birth date drops a node out of the competition integral", {
  build <- function(birth_date, times) {
    ctrl <- Control()
    ctrl$node_density_in_birth_date <- birth_date
    st <- FF16_Strategy()
    st$control <- ctrl
    s <- Species("FF16", "FF16_Env")(st)
    env <- Environment("FF16")
    env$set_fixed_environment(1.0, 200)
    for (t in times) {
      env$time <- t
      s$compute_rates(env, 1.0, 1.0)
      s$introduce_new_node()
    }
    s$heights <- c(9, 6, 4, 2)
    s
  }

  ok <- build(TRUE, c(0, 1, 2, 3))
  bad <- build(TRUE, c(0, 1, 1, 3))
  ## Same nodes, same heights, same carried densities: only the grid differs.
  expect_equal(bad$heights, ok$heights)
  expect_equal(bad$log_densities_state, ok$log_densities_state)

  ## The monotonicity check is what catches it, and it does.
  expect_strictly_increasing(ok$node_times)
  expect_failure(expect_strictly_increasing(bad$node_times))

  ## The integral, meanwhile, silently loses 15%+ of the competition -- and the
  ## axis assertion above cannot see that, because both the C++ and the
  ## by-hand trapezium walk the same defective grid.
  expect_gt(abs(bad$compute_competition(0) - ok$compute_competition(0)) /
              ok$compute_competition(0), 0.1)
  expect_integrates_over(bad, quadrature_grid(bad), 0)

  ## The height coordinate does not integrate over these times, so it is
  ## unaffected: same nodes, same answer.
  expect_identical(build(FALSE, c(0, 1, 1, 3))$compute_competition(0),
                   build(FALSE, c(0, 1, 2, 3))$compute_competition(0))
})

## Defect: the wrong density equation for the coordinate. A density in birth
## date changes by mortality alone -- nothing moves a cohort along the birth-date
## axis -- while a density in height also compresses as the spacing between
## neighbouring sizes changes. Keeping the compression term on the birth-date
## path, or dropping it from the height path, is an exact per-node identity to
## check rather than something to infer from a trend.
test_that("only the height coordinate's density rate carries a compression term", {
  for (x in c("FF16", "K93")) {
    nm <- Node(x, environment_type(x))(
      switch(x, FF16 = FF16_Strategy(), K93 = K93_Strategy()))$ode_names
    i_mort <- which(nm == "mortality")
    i_dens <- which(nm == "log_density")

    for (birth_date in c(FALSE, TRUE)) {
      scm <- short_run(x, birth_date)
      s <- scm$patch$species[[1]]
      env <- scm$patch$environment
      compression <- 0
      for (j in seq_len(min(s$size, 25L))) {
        nd <- s$node_at(j)
        rates <- nd$ode_rates
        grad <- nd$growth_rate_gradient(env)
        compression <- max(compression, abs(grad))
        expect_identical(rates[[i_dens]],
                         if (birth_date) -rates[[i_mort]]
                         else -rates[[i_mort]] - grad)
      }
      ## The term is there to be got wrong: ~5e-2 per year on FF16, so carrying
      ## it on the birth-date path would not be a rounding difference.
      expect_gt(compression, 1e-3)
    }
  }
})

## The height coordinate's boundary condition is N(H_0) = birth_rate * pr_estab
## / g(H_0) (eq-bc1); the birth-date one drops the division, because nu is a
## density per unit birth date and offspring arrive at a rate. Recording g(H_0)
## lets each be checked against the other's.
test_that("each coordinate's boundary condition is the one its integral needs", {
  for (birth_date in c(FALSE, TRUE)) {
    scm <- collected_run(birth_date)
    nd <- scm$patch$species[[1]]$new_node
    g0 <- nd$growth_rate_at_birth
    expect_gt(g0, 0)
    ## log(birth_rate * pr_estab), which is what the birth-date path carries
    ## directly. Compare in log space.
    undivided <- log(1.0 * nd$individual$establishment_probability(
      scm$patch$environment))
    expect_equal(nd$log_density + if (birth_date) 0 else log(g0), undivided,
                 tolerance = 1e-10)
    ## And the division is doing real work, so the two are not interchangeable.
    expect_gt(abs(log(g0)), 1)
  }
})

## The boundary node's birth date is the current time, and compute_rates()
## (which stamps it) runs after the set_ode_state() that rebuilds the
## environment. Reading the stamp during that rebuild would use the previous
## derivs call's time, making the last trapezium interval a function of the
## step size. The measured effect on offspring production is small (<1e-6 on
## FF16), so assert the invariant directly rather than through a result.
test_that("the boundary node's birth date tracks patch time", {
  for (x in c("FF16", "K93")) {
    p <- size_only_parameters(x)
    ctrl <- Control()
    ctrl$node_density_in_birth_date <- TRUE
    patch <- Patch(x, environment_type(x))(p, Environment(x), ctrl)

    for (t in c(0.0, 1.5, 7.25)) {
      patch$set_time(t)
      patch$compute_environment()
      expect_identical(patch$species[[1]]$new_node$introduction_time, t)
    }
  }
})

## A scheduled run cannot produce repeated introduction times; a patch seeded
## without per-node times can, and used to do so silently. Given what that costs
## the integral (above), Patch rejects it up front.
test_that("nodes sharing a birth date are rejected in birth-date coordinates", {
  x <- "FF16"; e <- "FF16_Env"
  p <- size_only_parameters(x)

  scm <- SCM(x, e)(p, Environment(x), Control())
  scm$collect <- TRUE
  scm$run()
  state <- export_patch_state(scm, step = max(2L, length(scm$history) %/% 2L))

  ## A faithful resume carries per-node times, so both coordinates accept it.
  p_ok <- set_initial_state(p, state)
  ctrl_bd <- Control()
  ctrl_bd$node_density_in_birth_date <- TRUE
  expect_no_error(Patch(x, e)(p_ok, Environment(x), ctrl_bd))

  ## Drop them and every node inherits the boundary node's birth date.
  p_bad <- p_ok
  p_bad$initial_node_times <- numeric(0)
  expect_error(Patch(x, e)(p_bad, Environment(x), ctrl_bd),
               "sharing an introduction time")

  ## The height coordinate never integrates over these, so it is unaffected.
  expect_no_error(Patch(x, e)(p_bad, Environment(x), Control()))
})

## The coordinate is stored per strategy but cannot differ between the species
## of one patch: Patch::add_strategies() overwrites every strategy's control
## with the patch's. Pin that, because compute_competition() sums the species'
## contributions into a single optical depth and would otherwise be adding a
## birth-date integral to a height integral. This is also why the multi-species
## convergence check below needs only enough refinement to see a ratio: the
## coordinate cannot vary between species, so what is left to check there is the
## shared profile, not the axis.
test_that("the patch control decides the coordinate for every species", {
  x <- "FF16"; e <- "FF16_Env"
  ctrl_bd <- Control()
  ctrl_bd$node_density_in_birth_date <- TRUE

  p0 <- scm_base_parameters(x)
  p <- add_strategies(p0, trait_matrix(c(0.08, 0.2), "lma"),
                      birth_rate = list(1.0, 1.0))
  ## Deliberately disagree with the patch control; it should be ignored.
  p$strategies[[1]]$control <- Control()

  patch <- Patch(x, e)(p, Environment(x), ctrl_bd)
  for (s in patch$species) {
    expect_true(s$density_in_birth_date)
  }

  ## And the other way round: a patch built with the default control carries
  ## every species in height, whatever the strategies say.
  p$strategies[[1]]$control <- ctrl_bd
  patch_h <- Patch(x, e)(p, Environment(x), Control())
  for (s in patch_h$species) {
    expect_false(s$density_in_birth_date)
  }
})

## ---------------------------------------------------------------------------
## Reporting boundary: whichever coordinate the solver carries internally, R is
## handed a density in *height*, so tidy_outputs.R's `density`,
## interpolate_to_heights() and the plots keep their meaning. The carried
## quantity is reported alongside rather than lost.
## ---------------------------------------------------------------------------

## Compared in log space: |log N_b - log N_h| is the log of the density ratio,
## so 1e-3 means the densities agree to ~0.1%. A *relative* comparison of the logs
## is meaningless here because log_density crosses zero.
log_density_gap <- function(refine) {
  p_of <- function() {
    p0 <- scm_base_parameters("FF16")
    p0$max_patch_lifetime <- 20
    add_strategies(p0, trait_matrix(0.08, "lma"), birth_rate = 1.0)
  }
  one <- function(birth_date) {
    ctrl <- Control()
    ctrl$node_density_in_birth_date <- birth_date
    scm <- SCM("FF16", "FF16_Env")(interleave_schedule(p_of(), refine),
                                   Environment("FF16"), ctrl)
    scm$run()
    scm$patch$species[[1]]
  }
  sa <- one(FALSE)
  sb <- one(TRUE)
  list(gap = abs(sb$log_densities - sa$log_densities),
       state_gap = abs(sb$log_densities_state - sa$log_densities))
}

test_that("log_densities is a height density in both coordinates", {
  g0 <- log_density_gap(0)
  g1 <- log_density_gap(1)

  ## The typical node agrees closely, and refining the schedule improves it at
  ## about the 2nd order of the differenced Jacobian (measured ratios ~3.9).
  expect_lt(median(g0$gap), 5e-3)
  expect_gt(median(g0$gap) / median(g1$gap), 2.5)
  expect_gt(quantile(g0$gap, 0.9) / quantile(g1$gap, 0.9), 2.5)

  ## Deliberately not asserted on max(): the *worst* node does not converge.
  ## |dh/dtau| is a ratio of two differences that both shrink as the schedule is
  ## refined, so cancellation error sets a floor, and refinement adds nodes in
  ## the near-empty tail (log density ~ -11) where that floor is worst. The
  ## reconstructed height density is sound in aggregate and for plotting; it
  ## should not be trusted node-by-node out in the tail.

  ## And the conversion is doing real work: unconverted, the carried quantity
  ## differs from the height density by ~1 in log space, not ~1e-3.
  expect_gt(median(g0$state_gap), 0.5)
})

test_that("log_densities_state is the quantity actually integrated", {
  a <- collected_run(FALSE)
  b <- collected_run(TRUE)
  sa <- a$patch$species[[1]]
  sb <- b$patch$species[[1]]

  ## On the height path there is nothing to convert.
  expect_equal(sa$log_densities_state, sa$log_densities)
  ## On the birth-date path it is the raw per-node state, which Node reports
  ## unconverted (a lone node cannot form dh/dtau).
  expect_equal(sb$log_densities_state,
               vapply(sb$nodes, function(nd) nd$log_density, numeric(1)))
  expect_false(isTRUE(all.equal(sb$log_densities_state, sb$log_densities)))

  ## And the two are related by the Jacobian, boundary node last.
  jac <- sb$height_jacobian
  expect_length(jac, sb$size + 1L)
  expect_equal(sb$log_densities,
               sb$log_densities_state - log(head(jac, sb$size)))
})

## The reported height density is not merely self-consistent with the Jacobian:
## integrated over *height* it reproduces the competition integral the solver
## actually took over *birth date*. That closes the loop between the two
## coordinates on a single state, and it is where a Jacobian wrong by any factor
## shows up -- note that the central difference makes the change of variables
## exact term-by-term for interior nodes, so this is only a real constraint
## because it also covers the boundary node's exact g(H_0).
test_that("the reported height density integrates to the solver's competition", {
  for (birth_date in c(FALSE, TRUE)) {
    s <- collected_run(birth_date)$patch$species[[1]]
    jac <- s$height_jacobian
    ## log_densities is already the height density; the boundary node is not in
    ## it, so convert that one here.
    n <- c(exp(s$log_densities),
           exp(s$new_node$log_density -
                 if (birth_date) log(jac[[s$size + 1L]]) else 0))
    crown <- c(vapply(s$nodes,
                      function(nd) nd$individual$compute_competition(0),
                      numeric(1)),
               s$new_node$individual$compute_competition(0))
    expect_equal(trapezium_by_hand(quadrature_grid(s, FALSE), n * crown),
                 s$compute_competition(0), tolerance = 1e-3)
  }
})

## The boundary node needs no differencing: dh/dtau = -g(H_0) at birth, and
## compute_initial_conditions() re-evaluates that every step, so the recorded
## rate is both exact and current for that node alone.
test_that("the boundary node's Jacobian is exact, not differenced", {
  b <- collected_run(TRUE)
  sb <- b$patch$species[[1]]
  n <- sb$size

  g0 <- sb$new_node$growth_rate_at_birth
  expect_gt(g0, 0)
  expect_identical(sb$height_jacobian[[n + 1L]], g0)

  ## And it is a real improvement over the one-sided difference it replaces:
  ## the boundary node sits a whole introduction interval from its neighbour.
  one_sided <- abs((sb$new_node$height - sb$heights[[n]]) /
                   (sb$new_node$introduction_time - sb$node_times[[n]]))
  expect_gt(abs(one_sided - g0) / g0, 0.1)

  ## Interior nodes are unaffected -- they never used the boundary node's own
  ## Jacobian, only its (h, tau) as the far point of a central difference.
  expect_true(all(is.finite(sb$height_jacobian[seq_len(n)])))
})

test_that("the collected state carries both densities", {
  a <- collected_run(FALSE)
  b <- collected_run(TRUE)

  rows_h <- rownames(a$patch$state$species[[1]])
  rows_b <- rownames(b$patch$state$species[[1]])
  expect_false("log_density_state" %in% rows_h)
  expect_true("log_density_state" %in% rows_b)

  ## The matrix that feeds tidy_patch() agrees with the accessor, so the
  ## `density = exp(log_density)` downstream is in height units.
  sb <- b$patch$species[[1]]
  st <- b$patch$state$species[[1]]
  expect_equal(unname(st["log_density", seq_len(sb$size)]), sb$log_densities)
  expect_equal(unname(st["log_density_state", seq_len(sb$size)]),
               sb$log_densities_state)
})

test_that("resume reads the raw state, so it is unaffected by the conversion", {
  x <- "FF16"; e <- "FF16_Env"
  ctrl <- Control()
  ctrl$node_density_in_birth_date <- TRUE

  scm <- collected_run(TRUE)
  state <- export_patch_state(scm, step = max(2L, length(scm$history) %/% 2L))
  p2 <- set_initial_state(scm$parameters, state)

  scm2 <- SCM(x, e)(p2, Environment(x), ctrl)
  scm2$collect <- TRUE
  scm2$run()

  ## The seeded patch reproduces the exported ODE state exactly -- which it
  ## could not if the reporting conversion had leaked into the resume path.
  expect_equal(scm2$history[[1]]$ode_state, state$ode_state)
})

## ---------------------------------------------------------------------------
## End-to-end convergence. Everything above is an invariant of a single patch
## state; these check that the invariants compose, through a full run and a
## shared competition profile, into two coordinates that agree.
##
## These are the only expensive tests in the file, and the cost is cubic in the
## node count -- halving the node spacing costs ~9x, so the ladder is dominated
## entirely by its last rung. Measured seconds for both coordinates at 1x / 2x /
## 4x the default schedule: FF16 one species 0.5 / 3.6 / 41; K93 one species
## 0.2 / 0.8 / 15; FF16 two species 1.7 / 21 / 175; K93 three species
## 0.7 / 9 / 86. Each case below therefore runs the shortest ladder that still
## shows its own signal, and nothing is refined past 2x.
## ---------------------------------------------------------------------------

ff16_two_species_parameters <- function() {
  p0 <- scm_base_parameters("FF16")
  ## The published setup runs to the default patch lifetime (105 y). Shortening
  ## it to 20 y takes 1x + 2x from 23 s to 6.5 s and leaves the convergence
  ## *rate* intact -- per-halving gap ratios 3.9 and 4.1, against 3.9 and 3.1 at
  ## the full lifetime. What it does not do is leave the gap itself unchanged:
  ## the gap is generated during recruitment but offspring production keeps
  ## accumulating afterwards, so a longer run dilutes it (species 2's gap at 2x
  ## is 7.6e-3 here against 2.3e-3 at the full lifetime). The ratio is what this
  ## case is asserted on; K93 below carries the tight absolute bound.
  p0$max_patch_lifetime <- 20
  add_strategies(p0, trait_matrix(c(0.0825, 0.2625), "lma"),
                 hyperpar = FF16_hyperpar,
                 birth_rate = list(11.99177, 16.51006))
}

k93_three_species_parameters <- function() {
  p0 <- scm_base_parameters("K93")
  p0$max_patch_lifetime <- 35.10667
  sp <- trait_matrix(c(0.042, 0.063, 0.052,
                       8.5e-3, 0.014, 0.015,
                       2.2e-4, 4.6e-4, 3e-4,
                       0.008, 0.008, 0.008,
                       1.8e-4, 4.4e-4, 5.1e-4,
                       1.4e-4, 2.5e-3, 8.8e-3,
                       0.044, 0.044, 0.044),
                     c("b_0", "b_1", "b_2", "c_0", "c_1", "d_0", "d_1"))
  add_strategies(p0, sp, birth_rate = c(20, 20, 20))
}

convergence_cases <- list(
  FF16 = list(model = "FF16",
              parameters = function() size_only_parameters("FF16")),
  K93 = list(model = "K93",
             parameters = function() size_only_parameters("K93")),
  FF16_two_species = list(model = "FF16",
                          parameters = ff16_two_species_parameters),
  K93_three_species = list(model = "K93",
                           parameters = k93_three_species_parameters))

## Offspring production in both coordinates at a given refinement, computed once
## per (case, refinement). Several tests below read the same runs -- the
## convergence trend, the claim about which coordinate is converged, and the
## marginal-species check -- so sharing them is what keeps the second and third
## of those free.
offspring_pair <- local({
  cache <- list()
  function(case, refine) {
    key <- paste(case, refine, sep = "/")
    if (is.null(cache[[key]])) {
      x <- convergence_cases[[case]]$model
      p <- interleave_schedule(convergence_cases[[case]]$parameters(), refine)
      cache[[key]] <<- list(height = offspring_in_coordinate(p, x, FALSE),
                            birth = offspring_in_coordinate(p, x, TRUE))
    }
    cache[[key]]
  }
})

coordinate_gap <- function(case, refine) {
  o <- offspring_pair(case, refine)
  abs(o$birth - o$height) / abs(o$height)
}

## The load-bearing claim: for FF16 and K93, growth is a function of size, so
## the compression term the height coordinate computes is correct and the two
## coordinates must converge to the same offspring production. A stale or
## duplicated abscissa makes the gap sit still under refinement instead of
## shrinking, which is what makes this a regression test rather than a
## restatement of the implementation.
##
## One halving each. The previous three-level form asserted
## gaps[1]/gaps[3] > 4, which is an average of only >2 per halving; asserting
## the per-halving ratio directly is the same claim or stronger, and the second
## halving cost 15 s on K93 and 41 s on FF16 against 1.0 s and 4.1 s for the
## first two levels together.
test_that("size-only strategies converge across density coordinates", {
  ## FF16: measured gaps 1.02e-2 -> 4.18e-3, ratio 2.44.
  ff16 <- lapply(0:1, function(r) coordinate_gap("FF16", r))
  expect_lt(ff16[[2]], ff16[[1]])
  expect_gt(ff16[[1]] / ff16[[2]], 2)
  expect_lt(ff16[[2]], 5e-3)

  ## K93: measured gaps 2.64e-3 -> 5.90e-4, ratio 4.48 -- the second order of
  ## the trapezium rule. (The next halving gives 1.43e-4, ratio 4.12, so the
  ## rate holds; it is not asserted because it costs 15 s to observe.) Loose
  ## bounds: the point is convergence, not the exact rate.
  k93 <- lapply(0:1, function(r) coordinate_gap("K93", r))
  expect_lt(k93[[2]], k93[[1]])
  expect_gt(k93[[1]] / k93[[2]], 3)
  ## And they really are converging to each other, not merely narrowing.
  expect_lt(k93[[2]], 2e-3)
})

## The birth-date answer is the accurate one: it is what both coordinates
## converge to, and it gets there on a far coarser schedule. This is the
## substantive claim of the change, so pin it rather than leaving it in prose.
##
## Stated locally -- how far each coordinate moves when the schedule is halved --
## rather than against a converged reference, which needed an 8x schedule and
## 170 s for one run. The local form reuses the runs the test above already did.
## It is also specific to FF16: K93's two coordinates straddle the limit and its
## height answer happens to start marginally closer, so the claim is not a
## general one and is not asserted there.
test_that("the birth-date coordinate converges on a coarser schedule", {
  coarse <- offspring_pair("FF16", 0)
  finer <- offspring_pair("FF16", 1)

  move <- function(f) abs(coarse[[f]] - finer[[f]]) / abs(finer[[f]])
  ## Measured: the birth-date answer moves 1.8e-4 under a halving, the height
  ## answer 6.1e-3 -- 35x further, and still climbing.
  expect_lt(move("birth"), 1e-3)
  expect_gt(move("height") / move("birth"), 10)
})

## Multiple species share one competition profile: Patch::compute_competition()
## sums their per-species integrals into a single optical depth, each taken over
## that species' own node list. So the coordinate change has to hold for each
## species *through* a profile the others also shape. These use the established
## multi-species setups from test-strategy-ff16.R and test-strategy-k93.R.
##
## One halving only. An integral over a wrong or duplicated axis gives a gap
## that barely moves under refinement, so one halving separates it from the ~4
## per halving a working pair of coordinates gives; and the axis itself is
## already pinned exactly, per species, by the invariants at the top of the file.
test_that("multi-species runs converge across density coordinates", {
  ## Measured per-halving ratios: FF16 3.9 and 4.1, K93 4.1, 4.1 and 4.1.
  for (case in c("FF16_two_species", "K93_three_species")) {
    coarse <- coordinate_gap(case, 0)
    finer <- coordinate_gap(case, 1)
    ## Every species converges, not just the aggregate.
    expect_true(all(finer < coarse))
    expect_true(all(coarse / finer > 2))
    expect_true(all(finer < 1e-2))
  }
  ## And K93's, at its full lifetime, converge to within the same bound the
  ## single-species cases meet.
  expect_true(all(coordinate_gap("K93_three_species", 1) < 2e-3))
})

## K93's three-species case largely excludes its first species -- its offspring
## production is ~90x below the dominant one. A marginal species is where a
## change to the shared profile could bite hardest, so check its *relative*
## error is no worse than the dominant species', rather than only that the
## totals agree. Reads the runs the test above already did.
test_that("a competitively excluded species is no worse off", {
  o <- offspring_pair("K93_three_species", 1)
  h <- o$height
  rel <- coordinate_gap("K93_three_species", 1)

  ## Species 1 really is marginal, so the premise of the test holds.
  expect_lt(h[[1]] / max(h), 0.05)
  ## And its relative error is of the same order as the rest.
  expect_lt(rel[[1]], 3 * max(rel[-1]))
})
