for (x in c("FF16", "K93")) {
  
  context(sprintf("Tidy-patch-%s", x))

  p0 <- scm_base_parameters(x)

  if(x == "FF16")
    p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0,  birth_rate_list=1.0)
  if (x == "K93")
    p1 <- expand_parameters(trait_matrix(0.059, "b_0"), p0,  birth_rate_list=1.0)

  env <- Environment(x)
  ctrl <- scm_base_control()

  expect_silent(
    results <- run_scm(p1, env, ctrl, collect = TRUE)
  )

  expect_contains(names(results), c("steps",  "n_spp", "species", "env", "offspring_production",   "net_reproduction_ratios", "p"))


  expect_is(results, "list")
  expect_is(results$steps, "data.frame")
  expect_is(results$species, "data.frame")
  expect_is(results$offspring_production, "numeric")
  expect_is(results$net_reproduction_ratios, "numeric")

  expect_is(results$env, "list")
  expect_is(results$env$light_availability, "data.frame")

  core_vars <- c("step", "time", "patch_density", "species", "node", "height", "mortality", "fecundity", "offspring_produced_survival_weighted", "log_density", "density")
  expect_contains(names(results$species), core_vars)
  
  # Test integration gives correct number, using artificial scenario with known answer

  ## set density = a*H^-b, height sequence 1:20, b = 3/2, a = starting density per unit area
  ## create stand of individuals with heights obtained and resulting density over time

  a <- 10
  b <- 3 / 2
  plant_stand <-
    tibble::tibble(
      time = 1,
      step = 2,
      patch_density = 1,
      species = 1,
      height = seq(20, 1, by = -0.01)
    ) %>%
    dplyr::mutate(
      node = seq_len(dplyr::n()) %>% rev(),
      density = a * height^(-b)
    )

  # integrate with plant inbuilt solver
  stand_integrate <- integrate_over_size_distribution(plant_stand)

  # Analyitcal soluntions
  # N: a/(1-b) H ^ (1-b)
  # H: a/(2-b) H ^ (2-b)

  N_analytical <- a / (1 - b) * 20^(1 - b) - a / (1 - b) * 1^(1 - b)
  N_trapezium <- -plant:::trapezium(plant_stand$height, plant_stand$density)
  N_plant <- stand_integrate$density

  expect_equal(N_analytical, N_trapezium, tolerance = 0.001)
  expect_equal(N_analytical, N_plant, tolerance = 0.001)
  
  H_analytical <- a / (2 - b) * 20^(2 - b) - a / (2 - b) * 1^(2 - b)
  H_trapezium <- -plant:::trapezium(plant_stand$height, plant_stand$height * plant_stand$density)
  H_plant <- stand_integrate$height

  expect_equal(H_analytical, H_trapezium, tolerance = 0.001)
  expect_equal(H_analytical, H_plant, tolerance = 0.001)
  
  Hav_analytical <- H_analytical / N_analytical
  Hav_trapezium <- H_trapezium / N_trapezium
  Hav_plant <- H_plant / N_plant

  expect_equal(Hav_analytical, Hav_trapezium, tolerance = 0.001)
  expect_equal(Hav_analytical, Hav_plant, tolerance = 0.001)


  # TODO: test tidy individual
  
  times <- c(1, 5, 10)
  expect_silent(
    tidy_species_new <- interpolate_to_times(results$species, times)
  )
  expect_true(all(names(tidy_species_new) %in% setdiff(names(results$species), c("step"))))
  # TODO: test interpolation actually works, also failures

  if(x == "FF16") {
    heights <- c(1, 5, 10)
    expect_silent(
      tidy_species_new <- interpolate_to_heights(results$species, heights)
    )
    expect_true(all(names(tidy_species_new) %in% setdiff(names(results$species), c("node"))))
    # TODO: test interpolation actually works, also failures
  }
}
