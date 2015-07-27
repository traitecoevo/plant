context("FFW16_PlantPlus")

test_that("Reference comparison", {
  cmp <- make_reference_plant()
  s <- FFW16_Strategy()
  p <- FFW16_PlantPlus(s)

  expect_that(p$strategy, is_identical_to(s))

  ## Set the height to something (here 10)
  h0 <- 10
  p$height <- h0

  vars <- p$internals

  expect_that(vars[["height"]],
              is_identical_to(h0))
  expect_that(vars[["area_leaf"]],
              equals(cmp$LeafArea(h0)))
  expect_that(vars[["mass_leaf"]],
              equals(cmp$LeafMass(cmp$traits$lma, cmp$LeafArea(h0))))
  expect_that(vars[["mass_sapwood"]],
              equals(cmp$SapwoodMass(cmp$traits$rho, cmp$LeafArea(h0), h0)))
  expect_that(vars[["mass_bark"]],
              equals(cmp$BarkMass(cmp$traits$rho, cmp$LeafArea(h0), h0)))
  expect_that(vars[["mass_root"]],
              equals(cmp$RootMass(cmp$LeafArea(h0))))
  expect_that(vars[["mass_live"]],
              equals(cmp$LiveMass(cmp$traits, cmp$LeafArea(h0))))
  expect_that(vars[["area_bark"]],
              equals(cmp$area_bark(h0)))
  expect_that(vars[["area_sapwood"]],
              equals(cmp$area_sapwood(h0)))
  expect_that(vars[["area_stem"]],
              equals(cmp$area_stem(h0)))

  expect_that(p$height,    is_identical_to(vars[["height"]]))
  expect_that(p$area_leaf, is_identical_to(vars[["area_leaf"]]))

  ## Heartwood function
  ## TODO: Check with Daniel here -- might need updating.
  ## Expect zero unless it has been set otherwise
  expect_that(vars[["area_heartwood"]],
              is_identical_to(0.0))
  HA0 <- 1e-3
  p$area_heartwood <- HA0
  ## TODO: This is due to issues with how we think about size; see notes
  ## in plant.cpp around set_height and set_mass_heartwood.  Note that
  ## when running as an ODE, this gives the *wrong answer*.
  h <- p$height
  p$height <- h + .1 # trick plant into recomputing all size variables
  p$height <- h
  vars <- p$internals
  expect_that(vars[["area_heartwood"]],
              is_identical_to(HA0))
  expect_that(vars[["area_stem"]],
              equals(cmp$area_stem(h0) + HA0))
  # set heartwood back at zero for subsequent tests
  p$area_heartwood <- 0

  env <- test_environment(h0)
  light_env <- attr(env, "light_env") # underlying function

  ## The R model computes A_lf * area_leaf * Y * c_bio, wheras we just
  ## compute A_lf; will have to correct some numbers.
  cmp_const <- s$Y * s$c_bio

  ## TODO: Check what growth variables look like before running
  ## through with environment.  Most are NA_real_, but some are not

  ## Compute the physiological variables and retreive them.
  p$compute_vars_phys(env)
  p$compute_vars_growth() # NOTE: Compute immediately *after* vars_phys
  vars <- p$internals

  ## 1. Assimilation:
  cmp_assimilation_plant <- cmp$assimilation.plant(h0, light_env)
  expect_that(vars[["assimilation"]],
              equals(cmp_assimilation_plant / cmp_const))

  ## 2. Respiration:
  cmp_respiration <- cmp$respiration.given.height(cmp$traits, h0)
  expect_that(vars[["respiration"]],
              equals(cmp_respiration / cmp_const))

  ## 3. Turnover:
  cmp_turnover <- cmp$turnover.given.height(cmp$traits, h0)
  expect_that(vars[["turnover"]],
              equals(cmp_turnover))

  ## 4. Net production:
  cmp_net_mass_production_dt <- cmp$net.production(cmp$traits, h0, light_env)
  expect_that(vars[["net_mass_production_dt"]],
              equals(cmp_net_mass_production_dt, tolerance=1e-7))

  ## 5. Reproduction fraction
  cmp_fraction_allocation_reproduction <-
    cmp$ReproductiveAllocation(cmp$traits$hmat,h0)
  expect_that(vars[["fraction_allocation_reproduction"]],
              equals(cmp_fraction_allocation_reproduction))

  ## 6. Fecundity rate
  cmp_fecundity_dt <- cmp$fecundity_dt(cmp$traits, h0, light_env)
  expect_that(vars[["fecundity_dt"]],
              equals(cmp_fecundity_dt, tolerance=1e-7))

  ## 8. Growth rate for height
  cmp_height_dt <- cmp$height.growth.dt(cmp$traits, h0, light_env)
  expect_that(vars[["height_dt"]],
              equals(cmp_height_dt, tolerance=1e-7))

  cmp_height_dt <-
    cmp$height.growth.dt.via.area.leaf(cmp$traits, h0, light_env)
  expect_that(vars[["height_dt"]],
              equals(cmp_height_dt, tolerance=1e-7))

  ## 9. Mortality rate
  cmp_mortality_dt <- cmp$mortality.dt(cmp$traits, h0, light_env)
  expect_that(vars[["mortality_dt"]],
              equals(cmp_mortality_dt))

  ## 10. Archietcural layout
  cmp_dheight_darea_leaf <- cmp$dHdA(cmp$LeafArea(h0))
  expect_that(vars[["dheight_darea_leaf"]],
              equals(cmp_dheight_darea_leaf))

  ## 11. Sapwood mass per leaf mass
  cmp_dmass_sapwood_darea_leaf<- cmp$sapwood.per.leaf.area(cmp$traits, h0)
  expect_that(vars[["dmass_sapwood_darea_leaf"]],
              equals(cmp_dmass_sapwood_darea_leaf))

  ## 12. Bark mass per leaf mass
  cmp_dmass_bark_darea_leaf <- cmp$bark.per.leaf.area(cmp$traits, h0)
  expect_that(vars[["dmass_bark_darea_leaf"]],
              equals(cmp_dmass_bark_darea_leaf))

  ## 12. Root mass per leaf mass
  cmp_dmass_root_darea_leaf <- cmp$root.per.leaf.area(cmp$traits, h0)
  expect_that(vars[["dmass_root_darea_leaf"]],
              equals(cmp_dmass_root_darea_leaf))

  ## 13. Leaf area growth rate
  cmp_area_leaf_dt <- cmp$area_leaf_dt(cmp$traits, h0, light_env)
  expect_that(vars[["area_leaf_dt"]],
              equals(cmp_area_leaf_dt, tolerance=1e-7))

  ## 14. sapwood area growth rate
  cmp_area_sapwood_dt <- cmp$area_sapwood_dt(cmp$traits, h0, light_env)
  expect_that(vars[["area_sapwood_dt"]],
              equals(cmp_area_sapwood_dt, tolerance=1e-7))

  ## 15. bark area growth rate
  cmp_area_bark_dt <- cmp$area_bark_dt(cmp$traits, h0, light_env)
  expect_that(vars[["area_bark_dt"]],
              equals(cmp_area_bark_dt, tolerance=1e-7))

  ## 16. heartwood area growth rate
  cmp_area_heartwood_dt <- cmp$area_heartwood_dt(cmp$traits, h0, light_env)
  expect_that(vars[["area_heartwood_dt"]],
              equals(cmp_area_heartwood_dt, tolerance=1e-7))

  ## 17. basal area growth rate
  cmp_area_stem_dt <- cmp$area_stem_dt(cmp$traits, h0, light_env)
  expect_that(vars[["area_stem_dt"]],
              equals(cmp_area_stem_dt, tolerance=1e-7))

  ## 18. change in basal diam per basal area
  cmp_ddiameter_stem_darea_stem <- cmp$ddiameter_stem_darea_stem(cmp$area_stem(h0))
  expect_that(vars[["ddiameter_stem_darea_stem"]],
              equals(cmp_ddiameter_stem_darea_stem, tolerance=1e-7))

  ## 18. basal diam growth rate
  cmp_diameter_stem_dt <- cmp$diameter_stem_dt(cmp$traits, h0, light_env)
  expect_that(vars[["diameter_stem_dt"]],
              equals(cmp_diameter_stem_dt, tolerance=1e-7))

  ## Check that height decomposition multiplies out to give right
  ## answer
  cmp <- prod(unlist(vars[c("dheight_darea_leaf",
                            "darea_leaf_dmass_live",
                            "fraction_allocation_growth")]),
              vars[[c("net_mass_production_dt")]])
  expect_that(vars[["height_dt"]],
              equals(cmp, tolerance=1e-7))
})

test_that("Seed bits", {
  ## Seed stuff:
  s <- FFW16_Strategy()
  seed <- FFW16_PlantPlus(s)
  env <- test_environment(10) # high enough
  light_env <- attr(env, "light_env") # underlying function

  ## Check that our root-finding succeeded and the leaf mass is correct:
  expect_that(seed$internals[["mass_live"]],
              equals(s$mass_seed, tolerance=1e-7))

  ## Check that the height at birth is correct.  These answers are
  ## actually quite different, which could come from the root finding?
  cmp <- make_reference_plant()
  expect_that(seed$height,
              equals(cmp$height.at.birth(cmp$traits), tolerance=1e-4))

  ## Then, check the germination probabilities in the current light
  ## environment:
  expect_that(seed$germination_probability(env),
              equals(cmp$germination.probability(cmp$traits, light_env),
                     tolerance=1e-5))
})

## TODO: Missing here: all the plant growing stuff.  Move that
## elsewhere.

test_that("Assimilation over distribution", {
  s1 <- FFW16_Strategy()
  p1 <- FFW16_PlantPlus(s1)

  c2 <- Control(plant_assimilation_over_distribution=TRUE)
  s2 <- FFW16_Strategy(control=c2)
  p2 <- FFW16_PlantPlus(s2)

  p1$height <- 10.0
  p2$height <- p1$height
  env <- test_environment(p1$height)

  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)
  p1_vars <- p1$internals
  p2_vars <- p2$internals

  ## Result is similar but not identical:
  expect_that(p2_vars, equals(p1_vars, tolerance=1e-7))
  expect_that(p2_vars, not(is_identical_to(p1_vars)))
})

test_that("Non-adaptive assimilation integration works", {
  c1 <- Control(plant_assimilation_adaptive=TRUE,
                plant_assimilation_over_distribution=TRUE)
  s1 <- FFW16_Strategy(control=c1)
  p1 <- FFW16_PlantPlus(s1)

  c2 <- Control(plant_assimilation_adaptive=FALSE,
                plant_assimilation_over_distribution=TRUE)
  s2 <- FFW16_Strategy(control=c2)
  p2 <- FFW16_PlantPlus(s2)

  p1$height <- 10.0
  p2$height <- p1$height
  env <- test_environment(p1$height)

  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)
  p1_vars <- p1$internals
  p2_vars <- p2$internals

  ## Result is similar but not identical:
  expect_that(p2_vars[["assimilation"]],
              equals(p1_vars[["assimilation"]], tolerance=1e-3))
  expect_that(p2_vars, not(is_identical_to(p1_vars)))
})

test_that("Ode interface", {
  p <- FFW16_PlantPlus(FFW16_Strategy())
  expect_that(p$ode_size, equals(5))
  expect_that(p$ode_state,
              equals(c(p$height, p$mortality, p$fecundity,
                       p$area_heartwood, p$mass_heartwood)))

  env <- test_environment(p$height * 10)
  p$compute_vars_phys(env)
  p$compute_vars_growth() # NOTE: Compute immediately *after* vars_phys
  expect_that(p$ode_state,
              equals(c(p$height, p$mortality, p$fecundity,
                       p$area_heartwood, p$mass_heartwood)))
  vars <- as.list(p$internals)
  expect_that(p$ode_rates,
              equals(c(vars$height_dt,
                       vars$mortality_dt,
                       vars$fecundity_dt,
                       vars$area_heartwood_dt,
                       vars$mass_heartwood_dt)))

  state_new <- c(p$height * 2, runif(p$ode_size - 1L))
  p$ode_state <- state_new
  expect_that(p$ode_state, is_identical_to(state_new))

  expect_that(p$ode_names,
              is_identical_to(c("height", "mortality", "fecundity",
                                "area_heartwood", "mass_heartwood")))
})

test_that("conversions", {
  cmp <- make_reference_plant()
  s <- FFW16_Strategy()
  p <- FFW16_Plant(s)
  p$height <- 10
  env <- test_environment(p$height)

  pp1 <- FFW16_plant_to_plant_plus(p, NULL)
  expect_that(pp1, is_a("PlantPlus<FFW16>"))
  expect_that(pp1$internals$assimilation, equals(NA_real_))
  expect_that(pp1$height, equals(p$height))

  pp2 <- FFW16_plant_to_plant_plus(p, env)
  expect_that(pp2$internals$assimilation, is_more_than(0))
  pp1$compute_vars_phys(env)
  pp1$compute_vars_growth()
  expect_that(pp1$internals, equals(pp2$internals))

  p2 <- pp1$to_plant()
  expect_that(p2, is_a("Plant<FFW16>"))
  expect_that(p2$internals$height_dt, equals(NA_real_))
  expect_that(p2$internals, equals(p$internals))
})
