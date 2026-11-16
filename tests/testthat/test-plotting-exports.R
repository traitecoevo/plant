## What a plotting package needs from a run and reads from plant rather than
## re-deriving: the crown, the size variable, K93's stem sizes and competition
## field, TF24's soil water potential, and the control a run ran under.

test_that("canopy_profile() is the Yokozawa profile at FF16's defaults", {
  u <- seq(0, 1, by = 0.05)
  crown <- canopy_profile(FF16_Model(), u)
  eta <- FF16_Model()$pars$eta
  expect_equal(crown$profile$u, u)
  expect_equal(crown$profile$Q, (1 - u^eta)^2, tolerance = 1e-14)
  expect_equal(crown$profile$q, 2 * eta * (1 - u^eta) * u^(eta - 1),
               tolerance = 1e-12)
  expect_equal(crown$eta_c, 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta))
  ## q is the leaf area's density in u: it integrates to the whole crown.
  total <- stats::integrate(function(x) canopy_profile(FF16_Model(), x)$profile$q,
                            0, 1)$value
  expect_equal(total, 1, tolerance = 1e-8)
})

test_that("canopy_profile() reads every model's own crown", {
  u <- c(0, 0.3, 0.7, 1)
  for (s in list(TF24_Model(), TF24f_Model(), TF24t_Model(), TF24w_Model(),
                 K93_Model())) {
    crown <- canopy_profile(s, u)
    eta <- s$pars$eta
    expect_equal(crown$profile$Q, (1 - u^eta)^2, tolerance = 1e-14)
    expect_equal(crown$eta_c, 1 - 2 / (1 + eta) + 1 / (1 + 2 * eta))
  }
  s <- FF16_Model()
  s$pars$eta <- 3
  expect_equal(canopy_profile(s, 0.5)$profile$Q, (1 - 0.5^3)^2)
  expect_error(canopy_profile(list(), 0.5), "Not a plant model")
  expect_error(canopy_profile(FF16_Model(), 1.5), "relative height")
  expect_error(canopy_profile(FF16_Model(), 0.5, control = list()), "Control")
})

test_that("a tidied run keeps the control that sets how its crowns cast shade", {
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- 2
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 5)
  ctrl <- control(shading_model = "flat-top-soft-box")
  out <- tidy(run(p, control = ctrl, record = "steps"))
  expect_identical(out$control, ctrl)
  u <- seq(0, 1, by = 0.1)
  s <- out$p$strategies[[1]]
  soft <- canopy_profile(s, u, control = out$control)
  ## The strategy alone carries no shading model, so it resolves to FF16's
  ## default, the smooth profile; the run's control gives the soft box.
  expect_equal(canopy_profile(s, u)$profile$Q, (1 - u^s$pars$eta)^2,
               tolerance = 1e-14)
  expect_false(isTRUE(all.equal(soft$profile$Q, (1 - u^s$pars$eta)^2)))
  expect_equal(soft$profile$Q[c(1, 11)], c(1, 0))
  expect_true(all(diff(soft$profile$Q) <= 0))
  ## q is the soft box's own density, -dQ/du, not the Yokozawa one.
  h <- 1e-6
  mid <- c(0.85, 0.9, 0.95)
  dQ <- (canopy_profile(s, mid + h, out$control)$profile$Q -
           canopy_profile(s, mid - h, out$control)$profile$Q) / (2 * h)
  expect_equal(canopy_profile(s, mid, out$control)$profile$q, -dQ, tolerance = 1e-6)
})

test_that("model_size() says what each model's size state measures", {
  expect_identical(model_size("FF16")[c("size", "size_units")],
                   list(size = "height", size_units = "m"))
  k93 <- model_size("K93")
  expect_identical(k93$size, "stem diameter")
  expect_identical(k93$size_units, "cm")
  expect_identical(k93$competitive_size, "basal area")
  ## A model that extends TF24 inherits its entries.
  expect_identical(model_size("TF24w")[-1], model_size("TF24")[-1])
  expect_identical(model_size(base_parameters("K93")), k93)
  expect_error(model_size("XX99"), "Unknown model")
  expect_error(model_size(list(time = 1)), "tidied deterministic run")
  ## The registry refuses a base model that does not say what its size is.
  base <- yaml::read_yaml(system.file("model-registry.yml", package = "plant"))
  base$models$K93$size_units <- NULL
  f <- withr::local_tempfile(fileext = ".yml")
  yaml::write_yaml(base, f)
  expect_error(suppressWarnings(model_registry(f)),
               "model 'K93' extends none but lacks size_units")
})

test_that("expand_state() pairs each species with its own strategy", {
  ## The species column is character, so from the tenth species on a split on
  ## it alone would order "10" before "2".
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- 2
  lma <- seq(0.06, 0.3, length.out = 11)
  p <- add_strategies(p, trait_matrix(lma, "lma"), birth_rate = rep(1, 11))
  sp <- expand_state(tidy(run(p, record = "steps")))$species
  sp <- sp[!is.na(sp$height), ]
  by_species <- tapply(sp$mass_leaf / sp$area_leaf, as.integer(sp$species), mean)
  expect_equal(unname(by_species), lma, tolerance = 1e-12)
})

test_that("K93's run gives stem sizes and records its competition field", {
  p <- testing_parameters("K93")
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(0.059, "b_0"), birth_rate = 20)
  out <- tidy(run(p, record = "steps"))
  sp <- expand_state(out)$species
  ok <- !is.na(sp$height)
  expect_equal(sp$diameter_stem[ok], sp$height[ok] / 100)
  expect_equal(sp$area_stem[ok], pi / 4 * (sp$height[ok] / 100)^2)
  light <- out$env$light_availability
  expect_true(all(c("height", "light_availability", "slope") %in% names(light)))
  last <- light[light$step == max(light$step), ]
  expect_gt(nrow(last), 1)
  expect_true(all(last$light_availability > 0 & last$light_availability <= 1))
})

test_that("TF24's record carries each layer's soil water potential", {
  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- 0.5
  p <- add_strategies(p, trait_matrix(0.07, "lma"))
  env <- Environment("TF24")
  env$set_soil_parameters(3, NULL, NULL, c(1e-3, 2e-3, 4e-3), NULL)
  env$set_soil_water_state(rep(0.3, 3))
  out <- tidy(run(p, env = env, record = "steps"))
  psi <- out$env$soil_psi
  moist <- out$env$soil_moist
  expect_equal(nrow(psi), nrow(moist))
  ## Layer 0 is what the environment's own conversion gives.
  first <- moist$step == 1
  expect_equal(psi$soil_psi[first][1], env$psi_from_soil_moist(moist$soil_moist[first][1]))
  ## At the start every layer holds the same moisture, so each layer's
  ## potential scales with its own a_psi.
  expect_equal(psi$soil_psi[first] / psi$soil_psi[first][1], c(1, 2, 4),
               tolerance = 1e-12)
})

test_that("recording the soil water potential counts no clamp", {
  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- 0.5
  p <- add_strategies(p, trait_matrix(0.07, "lma"))
  dry <- function() {
    env <- Environment("TF24")
    env$set_soil_water_state(rep(0.05, env$get_soil_number_of_depths()))
    env
  }
  clamps <- function(res) {
    stats::setNames(census_clamp_counts_tf24(res$solver)[[1]],
                    census_clamp_names_tf24())
  }
  recorded <- run(p, env = dry(), record = "steps")
  ## The run clamps the potential at its ceiling, and the record reads the
  ## capped value ...
  expect_gt(clamps(recorded)[["soil_potential_ceiling"]], 0)
  expect_equal(max(tidy(recorded)$env$soil_psi$soil_psi), 50)
  ## ... yet the run counts the same clamps as one that records nothing.
  expect_identical(clamps(recorded), clamps(run(p, env = dry())))
})
