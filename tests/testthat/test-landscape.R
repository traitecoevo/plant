# A landscape: patches, each a scenario, linked by flows and run through run().

# One TF24 patch under pulsed rain alone: a scenario on a soil column of the
# given layer widths. A thin surface layer holds little of a storm, so it
# sheds; the default column takes most of one in.
shallow_column <- c(0.05, 0.1, 0.15, 0.2, 0.25)
storm_times <- seq(0.25, 1.75, by = 0.25)
storm_depths <- c(0.05, 0.01, 0.08, 0.03, 0.06, 0.02, 0.07)

water_patch <- function(widths = NULL, lifetime = 2) {
  env <- Environment("TF24")
  if (!is.null(widths)) {
    env$set_soil_layer_widths(widths)
  }
  env$extrinsic_drivers_set_constant("rainfall", 0)
  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list()
  p <- add_strategies(p, trait_matrix(0.0825, "lma"))
  scenario(p, env, events = events(rainfall_pulse(storm_times, storm_depths)))
}

# A patch's water accounts at the end of its run, in metres.
water_accounts <- function(res) {
  flux <- res$solver$patch$environment$get_soil_water_state_cumulative_flux()
  list(rain = flux[[1]], infiltration = flux[[2]], shed = flux[[5]])
}

light_patch <- function(model, lifetime = 20) {
  p <- testing_parameters(model)
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list()
  traits <- switch(model, K93 = trait_matrix(0.059, "b_0"),
                   trait_matrix(0.0825, "lma"))
  scenario(add_strategies(p, traits))
}

test_that("unlinked patches run as each does alone, whatever their models", {
  land <- landscape(list(a = light_patch("FF16"), b = light_patch("K93")))
  expect_s3_class(land, "plant_landscape")
  expect_equal(nrow(land$flows), 0)

  res <- run(land)
  expect_s3_class(res, "plant_landscape_run")
  expect_named(res$patches, c("a", "b"))
  expect_equal(nrow(res$flows), 0)
  for (patch in names(land$patches)) {
    alone <- run(land$patches[[patch]])
    expect_identical(res$patches[[patch]]$offspring_production,
                     alone$offspring_production)
    expect_identical(res$patches[[patch]]$solver$patch$ode_state,
                     alone$solver$patch$ode_state)
  }

  ## An unnamed list is named by position.
  expect_named(landscape(unname(land$patches))$patches, c("1", "2"))
})

test_that("what a patch sheds reaches the patch below, and none is lost", {
  land <- landscape(list(ridge = water_patch(shallow_column),
                         gully = water_patch()),
                    resource_flow("ridge", "gully"))
  res <- run(land)
  ridge <- res$patches$ridge
  gully <- res$patches$gully

  ## The source runs as it does alone: nothing flows into it.
  expect_identical(ridge$offspring_production,
                   run(land$patches$ridge)$offspring_production)

  ## Every amount the ridge shed is a transfer, at the instant it was shed.
  log <- ridge$event_log
  shed <- vapply(log$applied, function(a) a[[2]], numeric(1))
  expect_gt(sum(shed > 0), 1)
  expect_equal(res$flows$time, log$time[shed > 0])
  expect_equal(res$flows$amount, shed[shed > 0])
  expect_equal(unique(res$flows$from), "ridge")
  expect_equal(unique(res$flows$to), "gully")

  ## The gully saw its own rain and the run-on, and took some of it in.
  alone <- run(land$patches$gully)
  expect_equal(water_accounts(gully)$rain,
               sum(storm_depths) + sum(res$flows$amount), tolerance = 1e-12)
  expect_gt(water_accounts(gully)$infiltration,
            water_accounts(alone)$infiltration + 0.01)

  ## Rain onto both patches went into one of them or left the gully.
  entered <- water_accounts(ridge)$infiltration + water_accounts(gully)$infiltration
  expect_equal(2 * sum(storm_depths), entered + water_accounts(gully)$shed,
               tolerance = 1e-12)

  ## The gully's scenario is what it ran: its storms and the run-on.
  expect_equal(length(gully$scenario$events$time),
               length(storm_times) + nrow(res$flows))
  ## The landscape is as given, so it runs again to the same result.
  expect_identical(res$landscape, land)
  expect_identical(run(res$landscape)$patches$gully$offspring_production,
                   gully$offspring_production)
})

test_that("a patch passes on what it sheds again, split by its flows' weights", {
  land <- landscape(list(gully = water_patch(),
                         mid = water_patch(shallow_column),
                         ridge = water_patch(shallow_column)),
                    resource_flow("ridge", c("mid", "gully"), weight = c(0.5, 0.25)),
                    resource_flow("mid", "gully"))
  ## Sources run first, whatever order the patches were given in.
  expect_equal(landscape_order(names(land$patches), land$flows),
               c("ridge", "mid", "gully"))
  res <- run(land)
  expect_named(res$patches, c("gully", "mid", "ridge"))

  accounts <- lapply(res$patches, water_accounts)
  passed <- tapply(res$flows$amount, list(res$flows$from, res$flows$to), sum)
  expect_equal(passed[["ridge", "mid"]], 0.5 * accounts$ridge$shed, tolerance = 1e-12)
  expect_equal(passed[["ridge", "gully"]], 0.25 * accounts$ridge$shed, tolerance = 1e-12)
  ## The middle patch sheds its own rain and the ridge's run-on alike.
  expect_equal(passed[["mid", "gully"]], accounts$mid$shed, tolerance = 1e-12)
  expect_gt(accounts$mid$shed, accounts$ridge$shed)

  ## A quarter of the ridge's run-off has no flow to follow and leaves, as
  ## does what the gully sheds.
  entered <- sum(vapply(accounts, function(a) a$infiltration, numeric(1)))
  left <- 0.25 * accounts$ridge$shed + accounts$gully$shed
  expect_equal(3 * sum(storm_depths), entered + left, tolerance = 1e-12)
})

test_that("a landscape is a value: saved and read back, it runs the same", {
  land <- landscape(list(ridge = water_patch(shallow_column),
                         gully = water_patch()),
                    resource_flow("ridge", "gully"))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(land, path)
  first <- run(land)
  again <- run(readRDS(path))
  expect_identical(again$flows, first$flows)
  expect_identical(again$patches$gully$offspring_production,
                   first$patches$gully$offspring_production)
})

test_that("a run's arguments reach every patch", {
  land <- landscape(list(ridge = water_patch(shallow_column),
                         gully = water_patch()),
                    resource_flow("ridge", "gully"))
  res <- run(land, record = "steps")
  out <- tidy(res)
  expect_named(out, c("ridge", "gully"))
  expect_true(all(c("species", "env") %in% names(out$gully)))

  ## Refinement runs a patch several times, and only its last run's shed
  ## amounts are passed on. Two default columns: the largest storms overfill
  ## those too.
  deep <- landscape(list(upper = water_patch(), lower = water_patch()),
                    resource_flow("upper", "lower"))
  refined <- run(deep, refine_schedule = TRUE)
  expect_gt(nrow(refined$flows), 0)
  expect_equal(sum(refined$flows$amount),
               water_accounts(refined$patches$upper)$shed, tolerance = 1e-12)
  expect_equal(refined$flows, run(deep)$flows)

  ## The stochastic solver's pulses shed as the deterministic one's do.
  set.seed(1)
  drawn <- run(land, method = "stochastic")
  expect_equal(drawn$method, "stochastic")
  expect_equal(sum(drawn$flows$amount),
               water_accounts(drawn$patches$ridge)$shed, tolerance = 1e-12)
  expect_gt(nrow(drawn$flows), 0)
})

test_that("a landscape that cannot be run is refused where it is made", {
  patches <- list(ridge = light_patch("FF16"), gully = light_patch("FF16"))

  expect_error(landscape(list()), "list of scenarios")
  expect_error(landscape(list(a = patches$ridge$p)), "list of scenarios")
  expect_error(landscape(stats::setNames(patches, c("a", "a"))),
               "name of its own")
  expect_error(landscape(list(a = light_patch("FF16", 20),
                              b = light_patch("FF16", 30))),
               "one span of time")
  expect_error(landscape(patches, data.frame(from = "ridge", to = "gully")),
               "resource_flow\\(\\) returns")

  expect_error(landscape(patches, resource_flow("ridge", "creek")),
               "\"creek\", which is not a patch")
  expect_error(landscape(patches, resource_flow("ridge", "gully"),
                         resource_flow("gully", "ridge")),
               "form a loop")
  expect_error(landscape(patches, resource_flow("ridge", "gully", weight = 0.6),
                         resource_flow("ridge", "gully", weight = 0.6)),
               "summing to more than 1")

  expect_error(resource_flow("ridge", "ridge"), "cannot flow into itself")
  expect_error(resource_flow("ridge", "gully", weight = 1.5), "between 0 and 1")
  expect_error(resource_flow("ridge", "gully", resource = 0), "counted from 1")
  expect_error(resource_flow(c("a", "b"), c("c", "d", "e")), "length 1 or 3")

  land <- landscape(patches)
  expect_error(run(land, until = 5), "cannot take `until` or `from`")
  expect_error(run(land, from = list()), "cannot take `until` or `from`")
  expect_error(run(land, times = 1), "does not take times")
})
