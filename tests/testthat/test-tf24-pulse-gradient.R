# The stand gradient across rain pulses. The sweep crosses a pulse by
# transposing the map the run applied at that instant. A pulse the layer takes
# whole adds a constant, so its rows are the identity's. A pulse that fills the
# layer and sheds the rest sets it to saturation whatever it held, so the
# layer's row is zero there, and infiltration and runoff move with the state it
# replaced.

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
  # At an ODE tolerance of 1e-6, where lma's column, which reaches the seed's
  # height, resolves; at the default its difference carries about 1e-5.
  run(p, pulse_environment(),
      Control(node_density_coordinate = "birth_date", ode_tol_rel = 1e-6,
              ode_tol_abs = 1e-6),
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

test_that("a pulse that fills the layer and sheds the rest is the pinned-step central difference", {
  # 0.5 m into a layer at 0.3: both pulses fill layer 1 and shed the rest. A
  # sweep that carried the filled layer's rows as the identity's misses this
  # difference by 6e-4 to 3.6e-3 at any ODE tolerance; held to 1e-6.
  depth <- 0.5
  solver <- pulse_run(pulse_stand(), depth)
  shed <- vapply(solver$event_log$applied, function(a) a[[2]], 0)
  expect_true(all(shed > 0))
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

test_that("the forward tangent walks the run's pulses", {
  # The ladder's trajectory reference steps the recording forward at an active
  # scalar, crossing each insertion by the same map the sweep transposes. Its
  # value is the run's census only if that map applies the pulses, and its
  # derivative is the sweep's only if the map carries a filled layer's rows.
  solver <- pulse_run(pulse_stand(), 0.5)
  traits <- census_trait_names(solver, "TF24")
  g <- stand_gradient(solver, traits = "1.lma", hyperpar = NULL)
  tangent <- ladder_trajectory_tangent_tf24(solver, as.numeric(traits == "1.lma"))
  metrics <- c(stand_census(solver), stand_fitness(solver))
  expect_identical(names(metrics), rownames(g$gradient))
  for (i in seq_along(metrics)) {
    m <- names(metrics)[[i]]
    expect_lt(abs(tangent$value[[i]] - metrics[[i]]), 1e-9 * abs(metrics[[i]]),
              label = paste("value of", m))
    expect_lt(abs(tangent$tangent[[i]] - g$gradient[m, "1.lma"]),
              1e-6 * abs(tangent$tangent[[i]]), label = paste("tangent of", m))
  }
})

test_that("the forward tangent walks a harvest between introductions", {
  # A harvest moves the vegetation, so the field the next step's first rates
  # are taken in is the one after it. At t = 1.3 no cohort is introduced, so the
  # instant is the harvest's alone.
  p <- pulse_stand()
  expect_false(any(abs(p$node_schedule_times[[1]] - 1.3) < 1e-12))
  solver <- run(p, pulse_environment(),
                Control(node_density_coordinate = "birth_date", ode_tol_rel = 1e-6,
                        ode_tol_abs = 1e-6),
                events = events(harvest(time = 1.3, fraction = 0.3)),
                record = "states")$solver
  traits <- census_trait_names(solver, "TF24")
  g <- stand_gradient(solver, traits = "1.lma", hyperpar = NULL)
  tangent <- ladder_trajectory_tangent_tf24(solver, as.numeric(traits == "1.lma"))
  metrics <- c(stand_census(solver), stand_fitness(solver))
  for (i in seq_along(metrics)) {
    m <- names(metrics)[[i]]
    expect_lt(abs(tangent$value[[i]] - metrics[[i]]), 1e-9 * abs(metrics[[i]]),
              label = paste("value of", m))
    expect_lt(abs(tangent$tangent[[i]] - g$gradient[m, "1.lma"]),
              1e-6 * abs(tangent$tangent[[i]]), label = paste("tangent of", m))
  }
})
