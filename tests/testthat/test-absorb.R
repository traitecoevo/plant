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
trapezium <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
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
    trapezium(c(sp$node_times, pt$time),
              c(exp(sp$log_densities_state), exp(sp$new_node$log_density))) +
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

test_that("an absorbed cohort is frozen in place", {
  ## In a continuous run a cohort leaves the integrated set one way only: by
  ## passing the absorbing height. Once absorbed it carries no density and its
  ## individual is no longer integrated, so its height never moves again.
  s <- new_scm(ff16(300), ctrl(absorb = 0.05))
  s$collect <- TRUE
  s$run()
  sp <- s$patch$species[[1]]
  absorbed <- which(sp$retired & sp$heights >= h_abs(0.05))
  expect_gt(length(absorbed), 10)
  expect_true(all(exp(sp$log_densities_state[absorbed]) == 0))
  for (j in absorbed[c(1, length(absorbed))]) {
    exists <- vapply(s$history, function(hh) hh$species[[1]]$size >= j, logical(1))
    h <- vapply(s$history[exists], function(hh) hh$species[[1]]$heights[[j]], numeric(1))
    r <- vapply(s$history[exists], function(hh) hh$species[[1]]$retired[[j]], logical(1))
    after <- which(r)
    expect_equal(h[after], rep(h[after[[1]]], length(after)))
    ## It was absorbed at or above the height, not before.
    expect_gte(h[after[[1]]], h_abs(0.05))
  }
})

test_that("cohorts an arrival window weights are never absorbed", {
  p <- window_regime(ff16(400), 50, 60)
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
  expect_error(stand_gradient(s, metrics = "leaf_area"), "not differentiable")
})
