context("TF24_plot_diagnostics")

test_that("TF24_plot_diagnostics assembles a TF24 stand diagnostic figure", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")

  # Minimal collected TF24 run carrying the soil-moisture environment series
  # and per-node leaf diagnostics that TF24_plot_diagnostics() consumes.
  max_patch_lifetime <- 2
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- max_patch_lifetime
  p1 <- expand_parameters(trait_matrix(0.07, "lma"), p0)

  env <- Environment("TF24")
  env$set_soil_number_of_depths(15)
  env$set_soil_water_state(rep(0.2, times = 15))
  x <- seq(0, max_patch_lifetime, length.out = 100)
  y <- 0.25 * sin(2 * pi * x) + 1
  env$extrinsic_drivers_set_variable("rainfall", x = x, y = y)
  ctrl <- Control()

  results <- run_scm(p1, env = env, ctrl = ctrl, collect = TRUE)

  p <- TF24_plot_diagnostics(results, x, y)

  expect_s3_class(p, "patchwork")
})
