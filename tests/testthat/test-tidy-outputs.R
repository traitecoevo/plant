for (x in c("FF16", "K93")) {
  
  context(sprintf("Tidy-patch-%s", x))

  p0 <- scm_base_parameters(x)

  if(x == "FF16")
    p1 <- expand_parameters(trait_matrix(0.08, "lma"), p0,  birth_rate_list=1.0)
  if (x == "K93")
    p1 <- expand_parameters(trait_matrix(0.059, "b_0"), p0,  birth_rate_list=1.0)

  env <- Environment(x)
  ctrl <- Control()

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


  times <- c(1, 5, 10)
  expect_silent(
    tidy_species_new <- interpolate_to_times(results$species, times)
  )
  expect_true(all(names(tidy_species_new) %in% setdiff(names(results$species), c("step"))))

  if(x == "FF16") {
    heights <- c(1, 5, 10)
    expect_silent(
      tidy_species_new <- interpolate_to_heights(results$species, heights)
    )
    expect_true(all(names(tidy_species_new) %in% setdiff(names(results$species), c("node"))))
  }
}

test_that("tidy_individual returns a tidy table of individual states over time", {
  ind <- FF16_Individual()
  env <- Environment("FF16")
  env$set_fixed_environment(1.0, 100)
  times <- seq(0, 50, length.out = 11)

  res <- grow_individual_to_time(ind, times, env)
  tidy <- tidy_individual(res)

  expect_s3_class(tidy, "tbl_df")
  expect_equal(nrow(tidy), length(times))
  expect_named(tidy, c("step", "time", "height", "mortality",
                       "fecundity", "area_heartwood", "mass_heartwood"))
  expect_equal(tidy$step, seq_along(times))
  expect_equal(tidy$time, times)
  # values match the raw solver output ...
  expect_equal(tidy$height, res$state[, "height"])
  # ... and the individual grows monotonically in full light
  expect_true(all(diff(tidy$height) > 0))
})

test_that("interpolate_to_times recovers known values and NAs out-of-range", {
  # two nodes with height exactly linear in time: node 1 -> 2t, node 2 -> 3t.
  # A natural spline through collinear points is the line itself, so the
  # interpolated values are exact.
  df <- tibble::tibble(
    species = 1,
    node    = rep(c(1, 2), each = 3),
    time    = rep(c(0, 1, 2), 2),
    height  = c(0, 2, 4, 0, 3, 6),
    density = 1
  )

  out <- interpolate_to_times(df, times = c(0.5, 1.5))
  expect_equal(out$height[out$node == 1], c(1.0, 3.0))
  expect_equal(out$height[out$node == 2], c(1.5, 4.5))
  expect_equal(out$time[out$node == 1], c(0.5, 1.5))

  # times outside the observed range return NA
  out_oor <- interpolate_to_times(df, times = c(-1, 3))
  expect_true(all(is.na(out_oor$height)))
})

test_that("interpolate_to_heights recovers known values, rebuilds density, appends largest node", {
  # leaf_area linear in height (10*h) and log_density linear (-(h-1)).
  df <- tibble::tibble(
    species     = 1,
    time        = 1,
    height      = c(1, 2, 3),
    log_density = c(0, -1, -2),
    leaf_area   = c(10, 20, 30)
  )

  # interpolated grid points, plus the largest individual (height 3) appended
  # rather than dropped (#352)
  out <- interpolate_to_heights(df, heights = c(1.5, 2.5))
  expect_equal(out$height, c(1.5, 2.5, 3))
  expect_equal(out$leaf_area, c(15, 25, 30))
  expect_equal(out$log_density, c(-0.5, -1.5, -2))
  # density is rebuilt as exp(log_density)
  expect_equal(out$density, exp(c(-0.5, -1.5, -2)))

  # grid points outside the observed range are dropped, but the largest
  # individual is still retained rather than discarded (#352)
  out_oor <- interpolate_to_heights(df, heights = c(0, 5))
  expect_equal(out_oor$height, 3)
  expect_equal(out_oor$leaf_area, 30)
})

test_that("interpolate_to_heights retains the largest size bracket on a coarse grid (#352)", {
  # reprex from the issue: tallest node (3.3) sits above the highest in-range
  # grid point (3), so the coarse grid would otherwise drop the largest class.
  data <- tibble::tibble(
    step        = 1,
    species     = 1,
    time        = 10,
    height      = c(0.1, 1.2, 2.3, 3.3),
    log_density = log(c(0.1, 0.1, 0.1, 0.2))
  )

  out <- interpolate_to_heights(data, heights = 1:4)
  # the actual largest individual is kept ...
  expect_true(3.3 %in% out$height)
  # ... the out-of-range grid point (4, above 3.3) is dropped ...
  expect_false(any(out$height == 4))
  # ... and no node is silently lost to NA
  expect_false(anyNA(out$height))
})
