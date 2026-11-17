# The stand gradient across rain pulses. The sweep crosses an event as if it
# left the state's rows unchanged, which is exact for a pulse the top layer takes
# whole. A pulse that fills the layer and sheds the rest sets it to saturation
# whatever it held, so its row there is zero and the sweep would carry one: that
# run is refused.

pulse_stand <- function() {
  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- 3
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
  full <- p$node_schedule_times[[1]]
  p$node_schedule_times <- list(full[round(seq(1, length(full), length.out = 12))])
  p
}

pulse_environment <- function() {
  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.3, env$get_soil_number_of_depths()))
  env$extrinsic_drivers_set_constant("rainfall", 1)
  env
}

pulse_run <- function(p, depth, record = TRUE) {
  run(p, pulse_environment(), Control(node_density_coordinate = "birth_date"),
      events = events(rainfall_pulse(time = c(1, 2), depth = depth)),
      record = if (record) "states" else character(0))$solver
}

test_that("a pulse the layer takes whole is the pinned-step central difference", {
  # 2 mm into a layer at 0.3: nothing is shed. Measured 1e-8 to 3e-8 against the
  # difference below at a relative step of 1e-5; held to 1e-6.
  depth <- 0.002
  solver <- pulse_run(pulse_stand(), depth)
  expect_true(all(vapply(solver$event_log$applied, function(a) a[[2]], 0) == 0))
  traits <- c("lma", "a_r1", "K_s")
  metrics <- c("leaf_area", "mass_above_ground")
  g <- stand_gradient(solver, metrics = metrics, traits = paste0("1.", traits),
                      hyperpar = NULL)
  pinned <- pin_steps(solver)
  census_of <- function(q) stand_census(pulse_run(q, depth, record = FALSE))[metrics]
  for (name in traits) {
    fd <- central_difference(pinned, 1, name, NULL, census_of, rel = 1e-5)
    for (m in metrics) {
      expect_lt(abs(g$gradient[m, paste0("1.", name)] - fd[[m]]), 1e-6 * abs(fd[[m]]),
                label = paste(m, name))
    }
  }
})

test_that("a pulse that fills the layer and sheds the rest is refused, by name", {
  solver <- pulse_run(pulse_stand(), 0.5)
  expect_true(any(vapply(solver$event_log$applied, function(a) a[[2]], 0) > 0))
  expect_error(stand_gradient(solver, metrics = "leaf_area", hyperpar = NULL),
               "filled layer 1 and shed .*has no row for a capped pulse")
})
