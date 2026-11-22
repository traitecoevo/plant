## A scenario keeps its environment as a recipe, so a saved scenario or result
## runs after readRDS(). Each test compares against the live environment, or
## the run in the session that made it, with identical(): the recipe adds no
## numerics, so nothing it rebuilds may differ in any bit.

## A TF24 environment set through every route a caller has: unequal layers,
## per-layer soil parameters, a starting water profile, scalar fields and a
## driver given as a series.
tf24_configured_environment <- function() {
  env <- Environment("TF24")
  set_tf24_soil(env, graded_widths(1.2, 4, 0.05), theta = c(0.30, 0.25, 0.2, 0.15),
                K_sat = c(150, 120, 90, 60))
  env$a_infil <- 0.8
  env$f_PAR <- 0.45
  env$extrinsic_drivers_set_constant("atm_vpd", 1.4)
  env$extrinsic_drivers_set_variable("rainfall", c(0, 1, 2, 3, 4, 5, 6, 7),
                                     c(0, 2, 2, 2, 0, 1, 3, 1))
  env
}

saved_and_read_back <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f))
  saveRDS(x, f)
  readRDS(f)
}

## Equal layers divide a stated total, which the running sum of their widths
## misses in the last bits for many columns, here 1.2 m in seven. The recipe
## keeps the total, so the rebuilt environment divides the same one again.
tf24_equal_layer_environment <- function() {
  env <- Environment("TF24")
  env$set_soil_layer_widths(c(0.6, 0.6))
  env$set_soil_number_of_depths(7)
  env
}

test_that("a TF24 environment's recipe rebuilds it exactly", {
  for (env in list(Environment("TF24"), tf24_configured_environment(),
                   tf24_equal_layer_environment())) {
    recipe <- environment_recipe(env)
    expect_s3_class(recipe, "plant_environment_recipe")
    rebuilt <- Environment(saved_and_read_back(recipe))
    expect_identical(environment_recipe(rebuilt), recipe)
    expect_identical(rebuilt$get_soil_depths(), env$get_soil_depths())
    expect_identical(rebuilt$get_soil_mid_depths(), env$get_soil_mid_depths())
    expect_identical(rebuilt$get_soil_water_state(), env$get_soil_water_state())
    expect_identical(rebuilt$depth, env$depth)
    rebuilt$set_soil_number_of_depths(3)
    env$set_soil_number_of_depths(3)
    expect_identical(rebuilt$get_soil_depths(), env$get_soil_depths())
    ## Between the knots too, where the slopes the series was given decide it.
    for (d in env$extrinsic_drivers_get_names()) {
      u <- seq(0, 7, length.out = 101)
      expect_identical(rebuilt$extrinsic_drivers_evaluate_range(d, u),
                       env$extrinsic_drivers_evaluate_range(d, u), info = d)
    }
  }
  ## The case is live: the total is not the bottom of the last layer.
  equal <- tf24_equal_layer_environment()
  expect_identical(equal$depth, 1.2)
  expect_false(identical(equal$get_soil_depths()[[7]], 1.2))
})

## Every run builds its environment from the recipe, so a setting the recipe
## leaves out never reaches a run, and nothing would say so. Every member R can
## reach on an environment is therefore named here as one of: kept by the
## recipe, reset by every run (Patch::reset clears the time, the light field
## and the soil water, and sets the shading model from the run's Control), or
## read-only. A member added to RcppR6_classes.yml fails this until it is put
## in one, and one the recipe keeps is added to recipe() and set_recipe(). The
## members are read off a live object, because R CMD check does not install
## the yml.
test_that("every environment member R can set is kept by the recipe or reset by a run", {
  members <- function(cl) {
    env <- switch(cl, FF16_Environment = FF16_Environment(),
                  K93_Environment = K93_Environment(),
                  TF24_Environment = TF24_Environment())
    setdiff(ls(env), c("clone", "initialize"))
  }
  reset <- c("time", "light_availability", "set_fixed_environment", "clear")
  read <- c("get_environment_at_height", "ode_size")
  kept <- list(
    FF16_Environment = character(0),
    K93_Environment = character(0),
    TF24_Environment = c(
      "soil_moist_sat", "K_sat", "water_mixing_rate", "a_psi", "n_psi",
      "a_infil", "b_infil", "f_PAR", "depth", "set_soil_number_of_depths",
      "set_soil_layer_widths", "set_soil_parameters", "set_soil_water_state",
      "extrinsic_drivers_set_constant", "extrinsic_drivers_set_variable",
      "recipe", "set_recipe"))
  reset_also <- list(
    FF16_Environment = c("set_shading_model", "compute_rates"),
    K93_Environment = character(0),
    TF24_Environment = c("add_resource_pulse", "add_water_pulse", "compute_rates"))
  read_also <- list(
    FF16_Environment = character(0),
    K93_Environment = character(0),
    TF24_Environment = c(
      "extrinsic_drivers_get_names", "extrinsic_drivers_evaluate",
      "extrinsic_drivers_evaluate_range", "psi_from_soil_moist",
      "soil_moist_from_psi", "get_soil_layer_widths", "get_soil_depths",
      "get_soil_mid_depths", "get_soil_number_of_depths",
      "get_soil_water_state", "get_soil_water_state_cumulative_flux",
      "get_SWdown", "get_atm_vpd", "get_ca", "get_leaf_temp",
      "get_atm_o2_kpa", "get_atm_kpa", "soil"))
  for (cl in names(kept)) {
    expect_setequal(members(cl), c(kept[[cl]], reset, reset_also[[cl]], read,
                                   read_also[[cl]]))
  }
  ## The fields the TF24 recipe keeps are its own names.
  recipe <- environment_recipe(Environment("TF24"))
  fields <- c("soil_moist_sat", "K_sat", "water_mixing_rate", "a_psi", "n_psi",
              "a_infil", "b_infil", "f_PAR", "depth")
  expect_true(all(fields %in% names(recipe)))
})

test_that("FF16 and K93 recipes hold their type alone", {
  for (model in c("FF16", "K93")) {
    recipe <- environment_recipe(Environment(model))
    expect_identical(unclass(recipe), list(type = environment_type(model)))
    expect_s3_class(Environment(recipe), class(Environment(model))[[1]])
  }
})

test_that("a scenario is a value: later changes to its env do not reach it", {
  p <- testing_parameters("TF24")
  env <- Environment("TF24")
  sc <- scenario(p, env)
  env$extrinsic_drivers_set_constant("rainfall", 0)
  expect_identical(sc$env$drivers$rainfall$value, 1)
  expect_identical(sc$env$soil_moist_sat, env$soil_moist_sat)
  expect_s3_class(scenario(p, sc$env)$env, "plant_environment_recipe")
})

test_that("a scenario refuses an environment of another model or a dead one", {
  expect_error(scenario(testing_parameters("FF16"), Environment("TF24")),
               "TF24_Env environment, and FF16 runs in a FF16_Env one")
  dead <- saved_and_read_back(Environment("TF24"))
  expect_error(scenario(testing_parameters("TF24"), dead),
               "read back with readRDS")
  expect_error(environment_recipe(list()), "must be an Environment")
})

test_that("a recipe whose depth disagrees with its widths is refused", {
  recipe <- environment_recipe(Environment("TF24"))
  recipe$depth <- recipe$depth + 0.01
  expect_error(Environment(recipe), "from the bottom of its last layer")
  recipe <- environment_recipe(Environment("FF16"))
  recipe$soil_moist_sat <- 0.4
  expect_error(Environment(recipe), "holds its type alone")
})

test_that("a scenario read back runs as the one that was saved", {
  p <- tf24_fitness_stand_parameters()
  sc <- scenario(p, tf24_configured_environment(),
                 Control(node_density_coordinate = "birth_date"))
  a <- run(sc)
  b <- run(saved_and_read_back(sc))
  expect_identical(b$solver$patch$ode_state, a$solver$patch$ode_state)
  expect_identical(b$offspring_production, a$offspring_production)

  p <- ff16_parameters(10)
  a <- run(p)
  b <- run(saved_and_read_back(a$scenario))
  expect_identical(b$solver$patch$ode_state, a$solver$patch$ode_state)
})

test_that("a result read back runs again to where it stood, and on", {
  p <- ff16_parameters(20)
  ev <- events(node_introductions(p), harvest(time = 8, fraction = 0.3))
  whole <- run(p, events = ev, record = "steps")
  part <- run(p, events = ev, until = 6, record = "steps")
  back <- saved_and_read_back(part)
  expect_error(tidy(back), "run\\(res, until = res\\$time\\)")

  again <- run(back, until = back$time)
  expect_identical(again$time, part$time)
  expect_identical(again$solver$patch$ode_state, part$solver$patch$ode_state)
  expect_identical(tidy(again)$species, tidy(part)$species)
  expect_error(run(back, until = 1), "past `until`")

  on <- run(back)
  expect_true(on$complete)
  expect_identical(on$solver$patch$ode_state, whole$solver$patch$ode_state)
  expect_identical(on$offspring_production, whole$offspring_production)
  expect_identical(on$event_log, whole$event_log)
})

test_that("a refined result read back takes the refined steps", {
  p <- ff16_parameters(10)
  a <- run(p, refine_schedule = TRUE)
  b <- run(saved_and_read_back(a), until = a$time)
  expect_identical(b$scenario$p, a$scenario$p)
  expect_identical(b$offspring_production, a$offspring_production)
})

test_that("an invaders' result read back is refused by name", {
  p <- ff16_parameters(10)
  res <- run(p)
  inv <- run_mutant(res, add_mutant(p, trait_matrix(0.09, "lma")))
  expect_error(run(saved_and_read_back(inv)), "run its resident's scenario again")
})

test_that("a recipe lacking a field keeps the fresh environment's value", {
  ## As a recipe saved before a field was added would.
  fresh <- environment_recipe(Environment("TF24"))
  configured <- environment_recipe(tf24_configured_environment())
  old <- configured
  old$water_mixing_rate <- NULL
  old$f_PAR <- NULL
  got <- environment_recipe(Environment(old))
  expect_identical(got$water_mixing_rate, fresh$water_mixing_rate)
  expect_identical(got$f_PAR, fresh$f_PAR)
  expect_identical(got[setdiff(names(got), c("water_mixing_rate", "f_PAR"))],
                   configured[setdiff(names(configured),
                                      c("water_mixing_rate", "f_PAR"))])
  ## A column given without its water starts at half saturation.
  no_water <- configured
  no_water$soil_water <- NULL
  expect_identical(Environment(no_water)$get_soil_water_state(),
                   rep(configured$soil_moist_sat * 0.5, 4))
})

test_that("a recipe saved with layer depths sets the column from its widths", {
  ## As a recipe saved when it carried the depths would, here with the depths
  ## equal layers were then given, built from the total.
  env <- tf24_equal_layer_environment()
  old <- environment_recipe(env)
  old$layer_depths <- (1:7) * (1.2 / 7)
  expect_false(identical(old$layer_depths, env$get_soil_depths()))
  rebuilt <- Environment(old)
  expect_identical(rebuilt$get_soil_depths(), env$get_soil_depths())
  expect_null(environment_recipe(rebuilt)$layer_depths)
})

test_that("a recipe with a field it does not know is refused by name", {
  recipe <- environment_recipe(Environment("TF24"))
  recipe$Ksat <- 100
  expect_error(Environment(recipe), "no field 'Ksat'")
  recipe <- environment_recipe(Environment("TF24"))
  recipe$layer_widths <- NULL
  expect_error(Environment(recipe), "without the layer_widths")
  recipe$depth <- NULL
  recipe$layer_depths <- (1:5) * 0.3
  expect_error(Environment(recipe), "without the layer_widths")
})

test_that("a layered recipe lacking a layer value takes the scalar field", {
  recipe <- environment_recipe(tf24_configured_environment())
  recipe$layered_soil_parameters$K_sat <- NULL
  got <- environment_recipe(Environment(recipe))$layered_soil_parameters
  expect_identical(got$K_sat, rep(recipe$K_sat, 4))
  expect_identical(got$n_psi, recipe$layered_soil_parameters$n_psi)
  expect_error(TF24_Environment()$set_recipe(list(1, 2)), "must be a named list")
})

test_that("a column of many equal layers rebuilds exactly", {
  env <- Environment("TF24")
  env$set_soil_number_of_depths(1000)
  rebuilt <- Environment(environment_recipe(env))
  expect_identical(rebuilt$get_soil_depths(), env$get_soil_depths())
})

test_that("a recipe serves where an environment is taken", {
  sc <- scenario(testing_parameters("TF24"), Environment("TF24"))
  expect_warning(rainfall_pulse(1, 0.01, env = sc$env), "continuous `rainfall`")
  p <- ff16_parameters(10)
  expect_identical(make_initial_state(p, heights = c(1, 2), densities = c(0.2, 0.1),
                                      env = scenario(p)$env),
                   make_initial_state(p, heights = c(1, 2), densities = c(0.2, 0.1)))
  ind <- FF16_Individual()
  expect_identical(run(ind, env = environment_recipe(FF16_Environment()),
                       times = c(0, 1, 2))$state,
                   run(ind, times = c(0, 1, 2))$state)
})

test_that("an invaders' result read back says how to rebuild it", {
  p <- ff16_parameters(10)
  res <- run(p)
  inv <- saved_and_read_back(run_mutant(res, add_mutant(p, trait_matrix(0.09, "lma"))))
  expect_error(run_mutant(inv, p), "resident's plant_run")
})
