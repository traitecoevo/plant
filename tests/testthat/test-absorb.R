## The absorbing top layer of a continuous patch (#519): a cohort that has grown
## to its strategy's asymptotic height is absorbed into one pinned top cohort per
## species, which carries its density, its share of the light field and uptake,
## and its future seed output. Without it a stand run on never equilibrates: the
## oldest cohorts edge toward their asymptote for as long as they carry density.

ff16 <- function(lifetime) {
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, lifetime, 1))
  p
}
ctrl <- function(absorb = 0) {
  Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb)
}
new_scm <- function(p, c) SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), c)
## Each node's density as it carries it: in birth date, birth rate times survival.
node_log_densities <- function(sp) {
  vapply(sp$nodes, function(n) n$log_density, numeric(1))
}
h_abs <- function(epsilon) {
  s <- FF16_Strategy()
  s$pars$hmat * (1 + log((1 - epsilon) / epsilon) / s$pars$a_f2)
}

test_that("the top layer is off by default and refused where it cannot work", {
  expect_equal(Control()$absorb_epsilon, 0)
  expect_false(run_scm(ff16(5))$patch$species[[1]]$top_layer$on)
  expect_error(run_scm(ff16(5), ctrl = Control(absorb_epsilon = 0.01)), "birth-date")
  expect_error(run_scm(ff16(5), ctrl = ctrl(absorb = 0.5)), "\\(0, 0.5\\)")
  ## K93 has no maturation height, so nothing says where its growth stops.
  p <- scm_base_parameters("K93")
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(0.059, "b_0"), birth_rate = 20)
  expect_error(run_scm(p, ctrl = ctrl(absorb = 0.01)), "asymptotic height")
})

test_that("absorbing a stand conserves its density to rounding at the instant", {
  ## A stand grown without absorption, so many cohorts stand above the height,
  ## loaded with absorption on: the next introduction absorbs all of them.
  burn <- new_scm(ff16(400), ctrl())
  burn$run_to(300)
  state <- export_patch_state(burn)
  state$ode_state <- c(state$ode_state, rep(0, 5 + 5))  # an empty top block
  patch <- Patch("FF16", "FF16_Env")(set_initial_state(ff16(400), state),
                                     Environment("FF16"), ctrl(absorb = 0.05))
  total <- function(pt) {
    sp <- pt$species[[1]]
    sum(sp$establishment_weights * exp(c(node_log_densities(sp),
                                         sp$new_node$log_density))) +
      sp$top_layer$density
  }
  sp <- patch$species[[1]]
  tall <- sum(sp$heights[-sp$size] >= h_abs(0.05))
  expect_gt(tall, 100)
  before <- total(patch)
  patch$introduce_new_node(1L, patch$time)
  sp <- patch$species[[1]]
  expect_equal(sum(sp$retired), tall)
  expect_gt(sp$top_layer$density, 0)
  expect_equal(sp$top_layer$height, h_abs(0.05))
  expect_gt(sp$top_layer$weight, 0)
  expect_equal(total(patch), before, tolerance = 1e-13)
})

test_that("a stand with a top layer reaches a fixed point close to the unbounded one", {
  free <- new_scm(ff16(600), ctrl())
  free$run()
  s <- new_scm(ff16(600), ctrl(absorb = 0.02))
  s$snapshot_interval <- 25
  s$run()
  sp <- s$patch$species[[1]]
  ## Bounded: the absorbed cohorts are frozen in place.
  expect_lt(sum(!sp$retired), 150)
  expect_gt(sp$top_layer$density, 0)
  expect_equal(sp$top_layer$height, h_abs(0.02))
  ## A fixed point: the top of the stand and its density stop moving.
  late <- tidy_snapshots(s)$species
  late <- late[late$time >= 300, ]
  expect_lt(diff(range(late$height_max)) / mean(late$height_max), 1e-7)
  expect_lt(diff(range(late$density)) / mean(late$density), 1e-8)
  ## And the stand it settles to is the unbounded stand's to a few parts per
  ## million, its seed output included.
  expect_equal(stand_census(s), stand_census(free), tolerance = 1e-5)
  expect_equal(s$net_reproduction_ratios, free$net_reproduction_ratios, tolerance = 1e-5)
})

test_that("a cohort is absorbed only at the absorbing height, and frozen there", {
  ## In a continuous run a cohort leaves the integrated set one way only: by
  ## passing the absorbing height. Every cohort that leaves it between two legs
  ## while still carrying density was absorbed, so it must stand at the height;
  ## and once out its height never moves again.
  s <- new_scm(ff16(300), ctrl(absorb = 0.05))
  s$collect <- TRUE
  s$run()
  hist <- lapply(s$history, function(hh) hh$species[[1]])
  checked <- 0
  for (r in seq_along(hist)[-1]) {
    before <- hist[[r - 1]]
    after <- hist[[r]]
    k <- seq_len(before$size)
    left <- which(!before$retired[k] & after$retired[k] &
                  exp(node_log_densities(before)[k]) > 1e-200)
    for (j in left) {
      expect_gte(after$heights[[j]], h_abs(0.05))
      later <- vapply(hist[r:length(hist)], function(h) h$heights[[j]], numeric(1))
      expect_equal(later, rep(later[[1]], length(later)))
      checked <- checked + 1
    }
  }
  expect_gt(checked, 10)
})

test_that("a run starts without the last run's top layer", {
  ## The layer belongs to a run: a second run on the same SCM -- refinement's,
  ## run_mutant()'s pinned resident rerun, a recording -- starts with none.
  s <- new_scm(ff16(150), ctrl(absorb = 0.05))
  s$snapshot_interval <- 50
  s$run()
  a <- s$net_reproduction_ratios
  y <- s$patch$ode_state
  expect_gt(s$patch$species[[1]]$top_layer$density, 0)
  s$run()
  expect_identical(s$net_reproduction_ratios, a)
  expect_identical(s$patch$ode_state, y)
  ## The rerun's first snapshot is the empty patch, with no layer pinning its top.
  first <- s$snapshots[[1]]$species[[1]]
  expect_equal(first$top_density, 0)
  expect_true(is.na(first$height_max) || first$height_max < 1)
})

test_that("an event reaches the trees the top layer holds", {
  ## A harvest of everything above 10 m takes the same fraction of the layer,
  ## pinned at the absorbing height, as of the canopy cohorts. Read just after
  ## it: a harvest of nothing a moment later puts a leg boundary there.
  p <- ff16(152)
  mark <- harvest(time = 150.51, fraction = 0)
  run_with <- function(ev) {
    s <- run_scm(p, ctrl = ctrl(absorb = 0.05), events = ev, run = FALSE)
    s$run_to(150)
    before <- s$patch$species[[1]]$top_layer
    s$run_to(150.51)
    list(before = before, after = s$patch$species[[1]]$top_layer)
  }
  none <- run_with(events(node_introductions(p), mark))
  cut <- run_with(events(node_introductions(p), mark,
                         harvest(time = 150.5, fraction = 0.9, size_min = 10)))
  expect_gt(none$before$density, 0)
  ## With the harvest the layer holds a tenth of what it holds without.
  expect_equal(cut$after$density / none$after$density, 0.1, tolerance = 1e-3)
  expect_equal(cut$after$weight / none$after$weight, 0.1, tolerance = 1e-3)
  ## Its seeds are past output and are not harvested.
  expect_equal(cut$after$seeds, none$after$seeds, tolerance = 1e-3)
})

test_that("a layer whose density is spent empties rather than lingering", {
  ## The layer's density decays linearly, and at the bottom of the subnormal
  ## range a decay step rounds back to the value it started from. Left alone, a
  ## layer whose inflow stopped would keep the canopy top at the absorbing
  ## height for good; below the smallest normal double it is emptied at the
  ## next introduction instant, and its seeds stay.
  s <- new_scm(ff16(160), ctrl(absorb = 0.05))
  s$run_to(150)
  state <- export_patch_state(s)
  n <- state$n[[1]]
  top <- 8 * n + 5 + 1:3
  seeds <- 8 * n + 5 + 4
  state$ode_state[top] <- 5e-322
  ## And no cohort left to absorb: every exported one dead, its hazard at the
  ## sentinel the birth rate of one puts its density at.
  state$ode_state[8 * (seq_len(n) - 1) + 2] <- 750
  q <- set_initial_state(ff16(160), state)
  r <- run_scm(q, ctrl = ctrl(absorb = 0.05), events = state$events, run = FALSE)
  r$snapshot_interval <- 1
  r$run_to(152)
  layer <- r$patch$species[[1]]$top_layer
  expect_equal(layer$density, 0)
  expect_equal(layer$weight, 0)
  expect_equal(layer$seeds, state$ode_state[[seeds]], tolerance = 1e-12)
  expect_lt(r$patch$species[[1]]$height_max, h_abs(0.05))
})

test_that("a frozen cohort's height does not move the field", {
  ## A frozen cohort carries no density, so where its height sits must not
  ## matter -- not even out of the live cohorts' order, where a younger cohort
  ## that overtook an older one and then froze can put it. The fast field build
  ## once stopped its crossing search at such a cohort.
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 10
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  p$node_schedule_times <- list(0:4)
  pa <- Patch("FF16", "FF16_Env")(p, Environment("FF16"), ctrl(absorb = 0.05))
  for (t in 0:3) {
    pa$set_time(t)
    pa$introduce_new_node(1, t)
  }
  ## A node's hazard sets its density, at a birth rate of one; each interval
  ## between introductions a year apart holds a year's establishment, the
  ## newest none yet.
  node <- function(h, hazard, interval = c(1, 0.5)) {
    c(h, hazard, 0, 0, 0, 0, interval)
  }
  heights <- c(0.5, 1.2, 2, 2.5, 3, 3.4, 4)
  field <- function(frozen_at) {
    pa$set_ode_state(c(node(4.39, 0.01), node(1.5, 0.01), node(frozen_at, 750),
                       node(1.07, 0.01, c(0, 0)), rep(0, 10)), 3)
    vapply(heights, pa$environment$get_environment_at_height, numeric(1))
  }
  out_of_order <- field(3.5)
  expect_identical(pa$species[[1]]$retired, c(FALSE, FALSE, TRUE, FALSE))
  expect_identical(out_of_order, field(1.3))
})

test_that("refinement, a bad epsilon and a bad absorbing height are refused", {
  expect_error(run_scm(ff16(20), ctrl = ctrl(absorb = 0.05), refine_schedule = TRUE),
               "top layer")
  expect_error(run_scm(ff16(5), ctrl = ctrl(absorb = -0.1)), "zero \\(off\\)")
  expect_error(run_scm(ff16(5), ctrl = ctrl(absorb = NaN)), "zero \\(off\\)")
  p <- ff16(5)
  p$strategies[[1]]$pars$a_f2 <- 0
  expect_error(run_scm(p, ctrl = ctrl(absorb = 0.05)), "a_f2")
  ## A negative a_f2 can still give a positive height, but allocation to
  ## reproduction then falls with height, so there is no asymptote to absorb at.
  p$strategies[[1]]$pars$a_f2 <- -50
  expect_error(run_scm(p, ctrl = ctrl(absorb = 0.05)), "a_f2")
})

test_that("a cohort born with no density is refused rather than scored zero", {
  ## A zero birth rate is the old idiom for a rare invader. With the layer on its
  ## cohorts would be frozen at birth and score nothing per seed.
  p <- ff16(5)
  p$strategies[[1]]$birth_rate_y <- 0
  expect_error(run_scm(p, ctrl = ctrl(absorb = 0.05)), "birth rate zero")
  ## Off, it runs as it always has.
  expect_true(is.finite(run_scm(p, ctrl = ctrl())$net_reproduction_ratios))
})

test_that("tidying a collected run says what the node rows leave out", {
  ## The layer is no node, so tidy_patch() and every integral over its rows
  ## miss what it carries.
  expect_warning(run_scm(ff16(150), ctrl = ctrl(absorb = 0.05), collect = TRUE),
                 "top layer")
  expect_no_warning(run_scm(ff16(20), ctrl = ctrl(), collect = TRUE))
})

## Guards for any change to the fitness integral (#658 routes it through a plain
## trapezium): the layer's seeds are scaled by the birth rate at the birth dates
## its cohorts carry, as a node's output is.
test_that("offspring production scales the layer's seeds by the birth rate", {
  p <- ff16(150)
  p$strategies[[1]]$birth_rate_y <- 2
  s <- run_scm(p, ctrl = ctrl(absorb = 0.05))
  layer <- s$patch$species[[1]]$top_layer
  expect_gt(layer$seeds, 0)
  expect_equal(layer$seeds_birth, 2 * layer$seeds, tolerance = 1e-12)
  expect_equal(s$offspring_production, 2 * s$net_reproduction_ratios,
               tolerance = 1e-12)
})

test_that("each species has its own layer, at its own height", {
  ## Two strategies with different maturation heights, absorbed at their own
  ## heights, each conserving its density as the layer takes it.
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 200
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.0825), "lma"), birth_rate = list(1, 1))
  p$strategies[[2]]$pars$hmat <- 10
  p$node_schedule_times <- rep(list(node_schedule_times_periodic(0, 200, 1)), 2)
  s <- new_scm(p, ctrl(absorb = 0.05))
  s$run()
  layers <- lapply(s$patch$species, function(sp) sp$top_layer)
  expect_true(all(vapply(layers, function(l) l$density > 0, logical(1))))
  expect_equal(layers[[1]]$height, h_abs(0.05))
  s2 <- p$strategies[[2]]$pars
  expect_equal(layers[[2]]$height, 10 * (1 + log(0.95 / 0.05) / s2$a_f2))
})

test_that("TF24 absorbs at its own asymptotic height", {
  ## TF24 declares an asymptotic height from the same reproductive allocation
  ## curve; a short stand with a low maturation height reaches it.
  p <- scm_base_parameters("TF24")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 25
  p <- add_strategies(p, trait_matrix(c(0.0825, 2), c("lma", "hmat")),
                      hyperpar = TF24_hyperpar, birth_rate = 1)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 25, 1))
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 1)
  s <- run_scm(p, env, Control(node_density_in_birth_date = TRUE, absorb_epsilon = 0.05,
                               ode_step_size_max = 0.05))
  layer <- s$patch$species[[1]]$top_layer
  pars <- p$strategies[[1]]$pars
  expect_gt(layer$density, 0)
  expect_equal(layer$height, pars$hmat * (1 + log(0.95 / 0.05) / pars$a_f2))
  expect_true(is.finite(s$net_reproduction_ratios))
})

test_that("cohorts an arrival window weights are never absorbed", {
  p <- fitness_window(ff16(400), 50, 60)
  s <- run_scm(p, ctrl = ctrl(absorb = 0.05))
  sp <- s$patch$species[[1]]
  in_window <- sp$node_times >= 50 & sp$node_times <= 60
  ## The window's cohorts grew past the height and kept their own nodes; others
  ## at the same height went into the layer.
  expect_true(all(sp$heights[in_window] >= h_abs(0.05)))
  expect_false(any(sp$retired[in_window]))
  expect_gt(sum(sp$retired[!in_window]), 0)
  ## Nothing the window weights reached the layer.
  expect_equal(sp$top_layer$weight, 0)
})

test_that("a window resumed from a burn-in scores its own cohorts, not the layer's", {
  ## The burn-in runs without a window, so every cohort it absorbs carries a
  ## fitness weight of one into the exported layer. Resumed into a window, the
  ## ratio is the mean over the window's cohorts, none of which is absorbed.
  burn <- new_scm(ff16(200), ctrl(absorb = 0.05))
  burn$run_to(150)
  expect_gt(burn$patch$species[[1]]$top_layer$weight, 0)
  state <- export_patch_state(burn)
  p2 <- fitness_window(set_initial_state(ff16(300), state), 150, 200)
  p2$node_schedule_times <- list(node_schedule_times_periodic(150, 300, 1))
  s <- run_scm(p2, ctrl = ctrl(absorb = 0.05))
  sp <- s$patch$species[[1]]
  out <- sp$net_reproduction_ratio_by_node * p2$strategies[[1]]$pars$S_D
  k <- which(sp$node_times >= 150 & sp$node_times <= 200)
  expect_gt(sp$top_layer$seeds, 0)
  expect_equal(s$net_reproduction_ratios,
               sum(establishment_weights_within(sp, k) * out[k]) / 50,
               tolerance = 1e-12)
})

test_that("the top layer survives an export and a resume", {
  straight <- new_scm(ff16(300), ctrl(absorb = 0.05))
  straight$run()
  first <- new_scm(ff16(300), ctrl(absorb = 0.05))
  first$run_to(200)
  state <- export_patch_state(first)
  expect_gt(first$patch$species[[1]]$top_layer$density, 0)
  p2 <- set_initial_state(ff16(300), state)
  p2$node_schedule_times <- list(node_schedule_times_periodic(200, 300, 1))
  resumed <- new_scm(p2, ctrl(absorb = 0.05))
  resumed$snapshot_interval <- 1
  resumed$run()
  ## The layer the resumed run starts from is the one exported.
  expect_equal(resumed$snapshots[[1]]$species[[1]]$top_density,
               first$patch$species[[1]]$top_layer$density)
  expect_equal(resumed$patch$species[[1]]$top_layer$density,
               straight$patch$species[[1]]$top_layer$density, tolerance = 1e-3)
  expect_equal(resumed$net_reproduction_ratios, straight$net_reproduction_ratios,
               tolerance = 1e-3)
})

test_that("a gradient refuses the top layer", {
  s <- run_scm(ff16(10), ctrl = ctrl(absorb = 0.05), record_trajectory = TRUE)
  expect_error(stand_gradient(s, metrics = "leaf_area", hyperpar = NULL), "not differentiable")
})
