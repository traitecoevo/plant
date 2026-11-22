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
  ## The flows table has the same columns whatever is carried, and none here.
  expect_named(land$flows, c("from", "to", "type", "out_of", "into", "weight"))

  res <- run(land)
  expect_s3_class(res, "plant_landscape_run")
  expect_named(res$patches, c("a", "b"))
  expect_equal(nrow(res$transfers), 0)
  expect_named(res$transfers,
               c("time", "from", "to", "type", "out_of", "into", "amount"))
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
                    flow("ridge", "gully", "resource"))
  res <- run(land)
  ridge <- res$patches$ridge
  gully <- res$patches$gully

  ## The source runs as it does alone: nothing flows into it.
  expect_identical(ridge$solver$patch$ode_state,
                   run(land$patches$ridge)$solver$patch$ode_state)

  ## Every amount the ridge shed is a transfer, at the instant it was shed.
  log <- ridge$event_log
  shed <- vapply(log$applied, function(a) a[[2]], numeric(1))
  expect_gt(sum(shed > 0), 1)
  expect_equal(res$transfers$time, log$time[shed > 0])
  expect_equal(res$transfers$amount, shed[shed > 0])
  expect_equal(unique(res$transfers$from), "ridge")
  expect_equal(unique(res$transfers$to), "gully")

  ## The gully saw its own rain and the run-on, and took some of it in.
  alone <- run(land$patches$gully)
  expect_equal(water_accounts(gully)$rain,
               sum(storm_depths) + sum(res$transfers$amount), tolerance = 1e-12)
  expect_gt(water_accounts(gully)$infiltration,
            water_accounts(alone)$infiltration + 0.01)

  ## Rain onto both patches went into one of them or left the gully.
  entered <- water_accounts(ridge)$infiltration + water_accounts(gully)$infiltration
  expect_equal(2 * sum(storm_depths), entered + water_accounts(gully)$shed,
               tolerance = 1e-12)

  ## The gully met each amount at the instant the ridge shed it, into its
  ## surface, after its own storm at that instant: its log holds the storm and
  ## then the run-on, and at every other storm the storm alone.
  met <- gully$event_log
  expect_equal(unique(met$target_index), 1)
  asked <- vapply(met$requested, function(a) a[[1]], numeric(1))
  expected <- lapply(seq_along(storm_times), function(k) {
    arrived <- res$transfers$amount[res$transfers$time == storm_times[[k]]]
    data.frame(time = storm_times[[k]], amount = c(storm_depths[[k]], arrived))
  })
  expected <- do.call(rbind, expected)
  expect_identical(met$time, expected$time)
  expect_equal(asked, expected$amount)

  ## The gully's scenario is what it ran: its storms and the run-on.
  expect_equal(length(gully$scenario$events$time),
               length(storm_times) + nrow(res$transfers))
  ## The landscape is as given, so it runs again to the same result.
  expect_identical(res$landscape, land)
  expect_identical(run(res$landscape)$patches$gully$solver$patch$ode_state,
                   gully$solver$patch$ode_state)
})

test_that("a flow leaves the pool it names and reaches the pool it names", {
  ## The ridge takes water into its second layer too, more than that layer
  ## holds, so two of its pools shed at one instant.
  ridge <- water_patch(shallow_column)
  ridge$events <- events(ridge$events, resource_pulse(c(0.5, 1.5), 2, 0.2))
  land <- landscape(list(ridge = ridge, gully = water_patch()),
                    flow("ridge", "gully", "resource", out_of = 2, into = 3,
                                  weight = 0.5))
  res <- run(land)

  log <- res$patches$ridge$event_log
  shed <- vapply(log$applied, function(a) a[[2]], numeric(1))
  deep <- log$target_index == 2
  expect_true(all(shed[deep] > 0))
  expect_gt(sum(shed[!deep] > 0), 0)
  ## Only what the second layer shed is passed, and half of it.
  expect_equal(res$transfers$time, c(0.5, 1.5))
  expect_equal(res$transfers$amount, 0.5 * shed[deep])
  expect_equal(res$transfers$type, c("resource", "resource"))
  expect_equal(res$transfers$out_of, c(2L, 2L))
  expect_equal(res$transfers$into, c(3L, 3L))

  ## It reached the gully's third layer, and nothing else did.
  met <- res$patches$gully$event_log
  expect_equal(met$time[met$target_index == 3], c(0.5, 1.5))
  expect_equal(vapply(met$requested[met$target_index == 3],
                      function(a) a[[1]], numeric(1)),
               res$transfers$amount)
  expect_equal(sum(met$target_index == 1), length(storm_times))
})

test_that("a patch passes on what it sheds again, split by its flows' weights", {
  land <- landscape(list(gully = water_patch(),
                         mid = water_patch(shallow_column),
                         ridge = water_patch(shallow_column)),
                    flow("ridge", c("mid", "gully"), "resource", weight = c(0.5, 0.25)),
                    flow("mid", "gully", "resource"))
  ## Sources run first, whatever order the patches were given in.
  expect_equal(landscape_order(names(land$patches), land$flows),
               c("ridge", "mid", "gully"))
  res <- run(land)
  expect_named(res$patches, c("gully", "mid", "ridge"))

  accounts <- lapply(res$patches, water_accounts)
  passed <- tapply(res$transfers$amount, list(res$transfers$from, res$transfers$to), sum)
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
                    flow("ridge", "gully", "resource"))
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(land, path)
  first <- run(land)
  again <- run(readRDS(path))
  expect_identical(again$transfers, first$transfers)
  expect_identical(again$patches$gully$solver$patch$ode_state,
                   first$patches$gully$solver$patch$ode_state)
})

test_that("a run's arguments reach every patch", {
  land <- landscape(list(ridge = water_patch(shallow_column),
                         gully = water_patch()),
                    flow("ridge", "gully", "resource"))
  res <- run(land, record = "steps")
  out <- tidy(res)
  expect_named(out, c("ridge", "gully"))
  expect_true(all(c("species", "env") %in% names(out$gully)))

  ## Refinement runs a patch several times, and only its last run's shed
  ## amounts are passed on. Two default columns, which the largest storms
  ## overfill too, on a schedule coarse enough to be refined.
  coarse <- water_patch()
  coarse$p$node_schedule_times <- list(seq(0, 2, length.out = 8))
  deep <- landscape(list(upper = coarse, lower = coarse),
                    flow("upper", "lower", "resource"))
  refined <- run(deep, refine_schedule = TRUE)
  upper <- refined$patches$upper
  expect_gt(length(upper$scenario$p$node_schedule_times[[1]]), 8)
  shed <- vapply(upper$event_log$applied, function(a) a[[2]], numeric(1))
  expect_gt(sum(shed > 0), 0)
  expect_equal(nrow(refined$transfers), sum(shed > 0))
  expect_equal(sum(refined$transfers$amount), water_accounts(upper)$shed,
               tolerance = 1e-12)

  ## The stochastic solver's pulses shed as the deterministic one's do.
  set.seed(1)
  drawn <- run(land, method = "stochastic")
  expect_equal(drawn$method, "stochastic")
  expect_equal(sum(drawn$transfers$amount),
               water_accounts(drawn$patches$ridge)$shed, tolerance = 1e-12)
  expect_gt(nrow(drawn$transfers), 0)
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
               "flow\\(\\) returns")

  expect_error(landscape(patches, flow("ridge", "creek", "resource")),
               "\"creek\", which is not a patch")
  expect_error(landscape(patches, flow("ridge", "gully", "resource"),
                         flow("gully", "ridge", "resource")),
               "form a loop")
  expect_error(landscape(patches, flow("ridge", "gully", "resource", weight = 0.6),
                         flow("ridge", "gully", "resource", weight = 0.6)),
               "summing to more than 1")

  expect_error(flow("ridge", "ridge", "resource"), "cannot flow into itself")
  expect_error(flow("ridge", "gully", "resource", weight = 1.5), "between 0 and 1")
  expect_error(flow("ridge", "gully", "resource", out_of = 0),
               "`out_of` and `into` are indices, counted from 1")
  expect_error(flow("ridge", "gully", "resource", into = Inf),
               "`out_of` and `into` are indices, counted from 1")
  expect_error(flow("ridge", "gully", "resource", out_of = 1e10),
               "`out_of` and `into` are indices, counted from 1")
  ## Shares out of different pools of one patch are not summed together.
  expect_no_error(landscape(patches,
                            flow("ridge", "gully", "resource", out_of = 1:2)))
  ## The kind is named, and a value no kind takes is refused, so a misspelt
  ## argument cannot pass as a column.
  expect_error(flow("ridge", "gully"), "`type` is the kind of flow: one of \"resource\"")
  expect_error(flow("ridge", "gully", "pollen"), "`type` is the kind of flow")
  expect_error(flow("ridge", "gully", "resource", wieght = 0.5),
               "A resource flow takes no further values; it was given wieght")
  expect_error(flow("ridge", "gully", "resource", 1, 1, 1, 7),
               "it was given an unnamed value")
  expect_error(flow(c("a", "b"), c("c", "d", "e"), "resource"), "length 1 or 3")

  expect_equal(nrow(flow(character(0), character(0), "resource")), 0)

  land <- landscape(patches)
  expect_error(run(land, until = 5), "cannot take `until` or `from`")
  expect_error(run(land, from = list()), "cannot take `until` or `from`")
  expect_error(run(land, times = 1), "does not take times")
})

test_that("a landscape changed after it was made is checked again when run", {
  land <- landscape(list(ridge = water_patch(shallow_column),
                         gully = water_patch()),
                    flow("ridge", "gully", "resource"))

  heavier <- land
  heavier$flows$weight <- 2
  expect_error(run(heavier), "between 0 and 1")

  ## A kind of flow plant does not have is refused by name.
  pollen <- land
  pollen$flows$type <- "pollen"
  expect_error(run(pollen), "\"pollen\" is not a kind of flow; the kinds are resource")

  astray <- land
  astray$flows$to <- "creek"
  expect_error(run(astray), "\"creek\", which is not a patch")

  longer <- land
  longer$patches$gully <- water_patch(lifetime = 3)
  expect_error(run(longer), "one span of time")

  ## A change that leaves it sound is run: half the run-off, half the run-on.
  half <- land
  half$flows$weight <- 0.5
  expect_equal(run(half)$transfers$amount, 0.5 * run(land)$transfers$amount)
})

test_that("a flow into a pool the receiver lacks names the patch it stopped at", {
  ## A light-only patch has no resource pools to take run-on.
  ridge <- water_patch(shallow_column)
  dry <- light_patch("FF16", lifetime = 2)
  land <- landscape(list(ridge = ridge, slope = dry),
                    flow("ridge", "slope", "resource"))
  expect_error(run(land),
               "Patch \"slope\", with what flowed into it: .*0 resources")

  ## A flow out of a pool the source lacks passes nothing.
  none <- run(landscape(list(a = dry, b = dry), flow("a", "b", "resource")))
  expect_equal(nrow(none$transfers), 0)
})

test_that("the run reads a kind of flow only through its three functions", {
  ## Two made-up kinds beside the resource kind. Each lets go fixed amounts
  ## whatever the patch did, and each leaves what arrived as a note of its own
  ## on the receiving scenario, so what reached which kind can be read back.
  note <- function(name) {
    function(scenario, arrivals) {
      scenario[[name]] <- arrivals
      scenario
    }
  }
  kinds <- c(flow_kinds(), list(
    token = list(
      check = function(flows) if (any(flows$into == 99L)) stop("no token 99"),
      leaving = function(result, flow) {
        data.frame(time = c(1, 2), amount = c(1, 1) * flow$out_of)
      },
      arriving = note("tokens")),
    tag = list(
      check = function(flows) NULL,
      leaving = function(result, flow) data.frame(time = 1.5, amount = 10),
      arriving = note("tags"))))
  local_mocked_bindings(flow_kinds = function() kinds)
  flow_of <- function(type, from, to, ...) flow(from, to, type, ...)

  patches <- list(a = light_patch("FF16", 3), b = light_patch("FF16", 3),
                  c = light_patch("FF16", 3))
  ## Weights are a kind's own to bound: these are no shares, and the two
  ## kinds out of `a` are not summed together.
  land <- landscape(patches,
                    flow_of("tag", "a", "b", weight = 3),
                    flow_of("token", "a", c("b", "c"), out_of = c(1, 2),
                            into = c(1, 2), weight = c(0.5, 2)),
                    flow_of("token", "b", "c", weight = 1))
  res <- run(land)

  ## Every transfer carries its own kind's amounts, in time order, those at
  ## one instant in the order their patches ran.
  expect_equal(res$transfers$time, c(1, 1, 1, 1.5, 2, 2, 2))
  expect_equal(res$transfers$type,
               c("token", "token", "token", "tag", "token", "token", "token"))
  expect_equal(res$transfers$from, c("a", "a", "b", "a", "a", "a", "b"))
  expect_equal(res$transfers$to, c("b", "c", "c", "b", "b", "c", "c"))
  expect_equal(res$transfers$out_of, c(1L, 2L, 1L, 1L, 1L, 2L, 1L))
  expect_equal(res$transfers$into, c(1L, 2L, 1L, 1L, 1L, 2L, 1L))
  expect_equal(res$transfers$amount, c(0.5, 4, 1, 30, 0.5, 4, 1))

  ## Each kind was handed its own arrivals and no other's.
  b <- res$patches$b$scenario
  expect_equal(b$tokens$type, c("token", "token"))
  expect_equal(b$tokens$amount, c(0.5, 0.5))
  expect_equal(b$tags$type, "tag")
  expect_equal(b$tags$amount, 30)
  c_ran <- res$patches$c$scenario
  expect_equal(c_ran$tokens$from, c("a", "b", "a", "b"))
  expect_null(c_ran$tags)
  expect_null(res$patches$a$scenario$tokens)

  ## A kind's check sees its own rows alone.
  expect_error(landscape(patches, flow_of("token", "a", "b", into = 99)),
               "no token 99")
  expect_no_error(landscape(patches, flow_of("tag", "a", "b", into = 99)))
  ## The run's own limit holds for every kind: nothing may feed back.
  expect_error(landscape(patches, flow_of("tag", "a", "a")), "form a loop")

  ## The type names the kind, as text.
  odd <- land
  odd$flows$type <- factor(odd$flows$type)
  expect_error(run(odd), "`type` is the name of its kind")

  ## What a kind reports as having left must be finite times and finite,
  ## non-negative amounts; anything else is refused under the patch's name.
  for (bad in list(data.frame(time = 1, amount = NaN),
                   data.frame(time = 1, amount = -1),
                   data.frame(time = NA_real_, amount = 1),
                   data.frame(time = 1))) {
    kinds$tag$leaving <- function(result, flow) bad
    expect_error(run(land), "Patch \"a\" reported what left it by a tag flow")
  }
  kinds$tag$leaving <- function(result, flow) stop("nothing to read")
  expect_error(run(land),
               "Patch \"a\", in what left it by a tag flow: nothing to read")
  kinds$tag$leaving <- function(result, flow) data.frame(time = 1.5, amount = 10)
  kinds$tag$arriving <- function(scenario, arrivals) stop("nowhere to put it")
  expect_error(run(land),
               "Patch \"b\", with what flowed into it: nowhere to put it")
})

test_that("a kind may give its rows columns of its own", {
  kinds <- flow_kinds()
  kinds$resource$columns <- "size_min"
  local_mocked_bindings(flow_kinds = function() kinds)
  banded <- flow("a", "b", "resource", size_min = 2)
  expect_error(flow("a", "b", "resource", size_max = 2),
               "A resource flow takes size_min; it was given size_max")
  flows <- bind_flows(list(empty_flows(), flow("b", "c", "resource"), banded))
  expect_s3_class(flows, "plant_flows")
  expect_named(flows, c(names(empty_flows()), "size_min"))
  expect_equal(flows$size_min, c(NA, 2))
  patches <- list(a = light_patch("FF16", 3), b = light_patch("FF16", 3),
                  c = light_patch("FF16", 3))
  expect_equal(nrow(run(landscape(patches, flow("b", "c", "resource"), banded))$transfers), 0)
})

test_that("each patch has an area, by default equal", {
  patches <- list(a = light_patch("FF16"), b = light_patch("FF16"),
                  c = light_patch("FF16"))
  expect_equal(landscape(patches)$area, c(a = 1, b = 1, c = 1))
  expect_equal(landscape(patches, area = c(3, 2, 1))$area, c(a = 3, b = 2, c = 1))
  ## Named for the patches, in any order.
  expect_equal(landscape(patches, area = c(c = 1, a = 3, b = 2))$area,
               c(a = 3, b = 2, c = 1))
  ## Flows and areas together.
  land <- landscape(patches, flow("a", "b", "resource"), area = c(3, 2, 1))
  expect_equal(nrow(land$flows), 1)

  expect_error(landscape(patches, area = c(1, 2)), "one positive number for each")
  expect_error(landscape(patches, area = c(1, 2, 0)), "one positive number for each")
  expect_error(landscape(patches, area = c(1, 2, NA)), "one positive number for each")
  expect_error(landscape(patches, area = c(a = 1, b = 2, d = 3)),
               "`area` is named for a, b, d, and the patches are a, b, c")
  ## Checked again when run, as the rest of a landscape is.
  land$area <- c(1, -1, 1)
  expect_error(run(land), "one positive number for each")
})

test_that("what leaves a larger patch reaches a smaller one deeper", {
  land <- landscape(list(ridge = water_patch(shallow_column),
                         gully = water_patch()),
                    flow("ridge", "gully", "resource"),
                    area = c(ridge = 3, gully = 1))
  res <- run(land)
  equal <- run(landscape(land$patches, flow("ridge", "gully", "resource")))

  ## The ridge runs as before; three times its area sheds three times the
  ## water, onto one unit of gully.
  expect_equal(res$transfers$time, equal$transfers$time)
  expect_equal(res$transfers$amount, 3 * equal$transfers$amount)

  ## Water over the whole landscape, each patch's depths times its area: rain
  ## onto both went into one of them or left the gully.
  ridge <- water_accounts(res$patches$ridge)
  gully <- water_accounts(res$patches$gully)
  expect_equal(gully$rain, sum(storm_depths) + sum(res$transfers$amount),
               tolerance = 1e-12)
  expect_equal((3 + 1) * sum(storm_depths),
               3 * ridge$infiltration + 1 * gully$infiltration + 1 * gully$shed,
               tolerance = 1e-12)
})

test_that("a landscape's fitness weights each patch by its share of the area", {
  thin <- light_patch("FF16", 30)
  thick <- thin
  thick$p <- add_strategies(thin$p, trait_matrix(0.12, "lma"),
                            keep_existing = FALSE)
  land <- landscape(list(a = thin, b = thick), area = c(3, 1))
  res <- run(land)
  each <- vapply(res$patches, function(r) r$offspring_production, numeric(1))
  expect_gt(abs(each[["a"]] / each[["b"]] - 1), 0.01)
  expect_equal(res$offspring_production, (3 * each[["a"]] + each[["b"]]) / 4)
  ratios <- vapply(res$patches, function(r) r$net_reproduction_ratios, numeric(1))
  expect_equal(res$net_reproduction_ratios,
               (3 * ratios[["a"]] + ratios[["b"]]) / 4)
  ## Equal patches give the plain mean.
  expect_equal(run(landscape(land$patches))$offspring_production, mean(each))

  ## Species are matched by position, so patches holding different numbers
  ## of them have no landscape-wide value; nor has a stochastic run.
  two <- thin
  two$p <- add_strategies(thin$p, trait_matrix(0.12, "lma"))
  mixed <- run(landscape(list(a = thin, b = two)))
  expect_null(mixed$offspring_production)
  expect_null(mixed$net_reproduction_ratios)
  set.seed(1)
  expect_null(run(land, method = "stochastic")$offspring_production)
})

test_that("invaders run in every patch of a landscape", {
  thin <- light_patch("FF16", 30)
  thick <- thin
  thick$p <- add_strategies(thin$p, trait_matrix(0.12, "lma"),
                            keep_existing = FALSE)
  land <- landscape(list(a = thin, b = thick), area = c(3, 1))
  copy_of <- function(sc, lma) add_mutant(sc$p, trait_matrix(lma, "lma"))

  ## A copy of each patch's own resident, replayed in it, scores as the
  ## resident does, patch by patch and over the landscape.
  res <- run(land)
  copies <- run_mutant(res, list(a = copy_of(thin, 0.0825), b = copy_of(thick, 0.12)))
  expect_s3_class(copies, "plant_landscape_run")
  expect_named(copies$patches, c("a", "b"))
  expect_false(is.null(copies$mutants))
  expect_equal(copies$offspring_production / res$offspring_production, 1,
               tolerance = 1e-8)
  ## Invaders are run against a resident's result, not against invaders'.
  expect_error(run_mutant(copies, copy_of(thin, 0.0825)),
               "must be a resident landscape's run")

  ## One Parameters is used in every patch: the same invader meets both
  ## residents, and its landscape-wide value is the weighted mean of the two.
  res <- run(land)
  one <- run_mutant(res, copy_of(thin, 0.1))
  each <- vapply(one$patches, function(r) r$offspring_production, numeric(1))
  expect_gt(abs(each[["a"]] / each[["b"]] - 1), 0.01)
  expect_equal(one$offspring_production, (3 * each[["a"]] + each[["b"]]) / 4)

  res <- run(land)
  expect_error(run_mutant(res, list(a = copy_of(thin, 0.1))),
               "one for each patch, named as the patches are \\(a, b\\)")
})

test_that("an invader in a patch meets what flowed into it", {
  land <- landscape(list(ridge = water_patch(shallow_column),
                         gully = water_patch()),
                    flow("ridge", "gully", "resource"))
  res <- run(land)
  resident <- vapply(res$patches, function(r) r$offspring_production, numeric(1))
  copies <- run_mutant(res, add_mutant(land$patches$gully$p,
                                       trait_matrix(0.0825, "lma")))
  copy <- vapply(copies$patches, function(r) r$offspring_production, numeric(1))
  expect_true(all(resident > 0))
  expect_lt(max(abs(copy / resident - 1)), 1e-6)
  ## The gully's resident did better with the run-on than alone, and so the
  ## copy must have met the run-on to match it.
  alone <- run(land$patches$gully)$offspring_production
  expect_gt(abs(resident[["gully"]] / alone - 1), 1e-3)
})
