context("Plant")

test_that("Reference comparison", {
  cmp <- make_reference_plant()
  s <- Strategy()
  p <- Plant(s)

  expect_that(p$strategy, is_identical_to(s))

  ## Set the height to something (here 10)
  h0 <- 10
  p$height <- h0

  p_size <- p$vars_size

  expect_that(p_size[["height"]],
              is_identical_to(h0))
  expect_that(p_size[["leaf_area"]],
              equals(cmp$LeafArea(h0)))
  expect_that(p_size[["leaf_mass"]],
              equals(cmp$LeafMass(cmp$traits$lma, cmp$LeafArea(h0))))
  expect_that(p_size[["sapwood_mass"]],
              equals(cmp$SapwoodMass(cmp$traits$rho, cmp$LeafArea(h0), h0)))
  expect_that(p_size[["bark_mass"]],
              equals(cmp$BarkMass(cmp$traits$rho, cmp$LeafArea(h0), h0)))
  expect_that(p_size[["root_mass"]],
              equals(cmp$RootMass(cmp$LeafArea(h0))))
  expect_that(p_size[["live_mass"]],
              equals(cmp$LiveMass(cmp$traits, cmp$LeafArea(h0))))
  expect_that(p_size[["bark_area"]],
              equals(cmp$bark_area(h0)))
  expect_that(p_size[["sapwood_area"]],
              equals(cmp$sapwood_area(h0)))
  expect_that(p_size[["basal_area"]],
              equals(cmp$basal_area(h0)))

  expect_that(p$height,    is_identical_to(p_size[["height"]]))
  expect_that(p$leaf_area, is_identical_to(p_size[["leaf_area"]]))

  ## Heartwood function
  ## TODO: Check with Daniel here -- might need updating.
  ## Expect zero unless it has been set otherwise
  expect_that(p_size[["heartwood_area"]],
              is_identical_to(0.0))
  HA0 <- 1e-3
  p$heartwood_area <- HA0
  ## TODO: This is due to issues with how we think about size; see notes
  ## in plant.cpp around set_height and set_heartwood_mass.  Note that
  ## when running as an ODE, this gives the *wrong answer*.
  h <- p$height
  p$height <- h + .1 # trick plant into recomputing all size variables
  p$height <- h
  p_size <- p$vars_size
  expect_that(p_size[["heartwood_area"]],
              is_identical_to(HA0))
  expect_that(p_size[["basal_area"]],
              equals(cmp$basal_area(h0) + HA0))
  # set heartwood back at zero for subsequent tests
  p$heartwood_area <- 0


  env <- test_environment(h0)
  light_env <- attr(env, "light_env") # underlying function

  ## The R model computes A_lf * leaf_area * Y * c_bio, wheras we just
  ## compute A_lf; will have to correct some numbers.
  cmp_const <- s$Y * s$c_bio

  ## TODO: Check what p$vars_growth looks like before running through
  ## with environment.  Most are NA_real_, but some are not (some are
  ## computed on exit, but that's an implementation detail).

  ## Compute the physiological variables and retreive them.
  p$compute_vars_phys(env)
  p_phys <- p$vars_phys
  p_growth <- p$vars_growth

  ## 1. Assimilation:
  cmp_assimilation_plant <- cmp$assimilation.plant(h0, light_env)
  expect_that(p_phys[["assimilation"]],
              equals(cmp_assimilation_plant / cmp_const))

  ## 2. Respiration:
  cmp_respiration <- cmp$respiration.given.height(cmp$traits, h0)
  expect_that(p_phys[["respiration"]],
              equals(cmp_respiration / cmp_const))

  ## 3. Turnover:
  cmp_turnover <- cmp$turnover.given.height(cmp$traits, h0)
  expect_that(p_phys[["turnover"]],
              equals(cmp_turnover))

  ## 4. Net production:
  cmp_net_production <- cmp$net.production(cmp$traits, h0, light_env)
  expect_that(p_phys[["net_production"]],
              equals(cmp_net_production, tolerance=1e-7))

  ## 5. Reproduction fraction
  cmp_reproduction_fraction <-
    cmp$ReproductiveAllocation(cmp$traits$hmat,h0)
  expect_that(p_phys[["reproduction_fraction"]],
              equals(cmp_reproduction_fraction))

  ## 6. Fecundity rate
  cmp_fecundity_rate <- cmp$fecundity.rate(cmp$traits, h0, light_env)
  expect_that(p_phys[["fecundity_rate"]],
              equals(cmp_fecundity_rate, tolerance=1e-7))

  ## 7. Fraction of whole plant (mass) growth that is leaf.
  cmp_leaf_fraction <- cmp$leaf.fraction(cmp$traits, h0)
  expect_that(p_phys[["leaf_fraction"]],
              equals(cmp_leaf_fraction))

  ## 8. Growth rate for height
  cmp_height_growth_rate <- cmp$height.growth.rate(cmp$traits, h0, light_env)
  expect_that(p_phys[["height_growth_rate"]],
              equals(cmp_height_growth_rate, tolerance=1e-7))

  cmp_height_growth_rate <-
    cmp$height.growth.rate.via.mass.leaf(cmp$traits, h0, light_env)
  expect_that(p_phys[["height_growth_rate"]],
              equals(cmp_height_growth_rate, tolerance=1e-7))

  ## 9. Mortality rate
  cmp_mortality_rate <- cmp$mortality.rate(cmp$traits, h0, light_env)
  expect_that(p_phys[["mortality_rate"]],
              equals(cmp_mortality_rate))

  ## 10. Archietcural layout
  cmp_dheight_dleaf_area <- cmp$dHdA(cmp$LeafArea(h0))
  expect_that(p_growth[["dheight_dleaf_area"]],
              equals(cmp_dheight_dleaf_area))

  ## 11. Sapwood mass per leaf mass
  cmp_dsapwood_mass_dleaf_mass <- cmp$sapwood.per.leaf.mass(cmp$traits, h0)
  expect_that(p_growth[["dsapwood_mass_dleaf_mass"]],
              equals(cmp_dsapwood_mass_dleaf_mass))

  ## 12. Bark mass per leaf mass
  cmp_dbark_mass_dleaf_mass <- cmp$bark.per.leaf.mass(cmp$traits, h0)
  expect_that(p_growth[["dbark_mass_dleaf_mass"]],
              equals(cmp_dbark_mass_dleaf_mass))

  ## 12. Root mass per leaf mass
  cmp_droot_mass_dleaf_mass <- cmp$root.per.leaf.mass(cmp$traits, h0)
  expect_that(p_growth[["droot_mass_dleaf_mass"]],
              equals(cmp_droot_mass_dleaf_mass))

  ## 13. Leaf area growth rate
  cmp_dleaf_area_dt <- cmp$dleaf_area_dt(cmp$traits, h0, light_env)
  expect_that(p_growth[["dleaf_area_dt"]],
              equals(cmp_dleaf_area_dt, tolerance=1e-7))

  ## 14. sapwood area growth rate
  cmp_dsapwood_area_dt <- cmp$dsapwood_area_dt(cmp$traits, h0, light_env)
  expect_that(p_growth[["dsapwood_area_dt"]],
              equals(cmp_dsapwood_area_dt, tolerance=1e-7))

  ## 15. bark area growth rate
  cmp_dbark_area_dt <- cmp$dbark_area_dt(cmp$traits, h0, light_env)
  expect_that(p_growth[["dbark_area_dt"]],
              equals(cmp_dbark_area_dt, tolerance=1e-7))

  ## 16. heartwood area growth rate
  cmp_dheartwood_area_dt <- cmp$dheartwood_area_dt(cmp$traits, h0, light_env)
  expect_that(p_growth[["dheartwood_area_dt"]],
              equals(cmp_dheartwood_area_dt, tolerance=1e-7))

  ## 17. basal area growth rate
  cmp_dbasal_area_dt <- cmp$dbasal_area_dt(cmp$traits, h0, light_env)
  expect_that(p_growth[["dbasal_area_dt"]],
              equals(cmp_dbasal_area_dt, tolerance=1e-7))

  ## 18. change in basal diam per basal area
  cmp_dbasal_diam_dbasal_area <- cmp$dbasal_diam_dbasal_area(cmp$basal_area(h0))
  expect_that(p_growth[["dbasal_diam_dbasal_area"]],
              equals(cmp_dbasal_diam_dbasal_area, tolerance=1e-7))

  ## 18. basal diam growth rate
  cmp_dbasal_diam_dt <- cmp$dbasal_diam_dt(cmp$traits, h0, light_env)
  expect_that(p_growth[["dbasal_diam_dt"]],
              equals(cmp_dbasal_diam_dt, tolerance=1e-7))

  ## Check that height decomposition multiplies out to give right answer
  expect_that(p_phys[["height_growth_rate"]],
              equals(
                prod(p_growth[c("dheight_dleaf_area","dleaf_area_dleaf_mass","growth_fraction")], p_phys[c("leaf_fraction","net_production")]), tolerance=1e-7))
})

test_that("Seed bits", {
  ## Seed stuff:
  s <- Strategy()
  seed <- Plant(s)
  env <- test_environment(10) # high enough
  light_env <- attr(env, "light_env") # underlying function

  ## Check that our root-finding succeeded and the leaf mass is correct:
  expect_that(seed$vars_size[["live_mass"]],
              equals(s$s, tolerance=1e-7))

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
  s1 <- Strategy()
  p1 <- Plant(s1)

  c2 <- Control(plant_assimilation_over_distribution=TRUE)
  s2 <- Strategy(control=c2)
  p2 <- Plant(s2)

  p1$height <- 10.0
  p2$height <- p1$height
  env <- test_environment(p1$height)

  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)
  p1_phys <- p1$vars_phys
  p2_phys <- p2$vars_phys

  ## Result is similar but not identical:
  expect_that(p2_phys, equals(p1_phys))
  expect_that(p2_phys, not(is_identical_to(p1_phys)))
})

test_that("Non-adaptive assimilation integration works", {
  c1 <- Control(plant_assimilation_adaptive=TRUE,
                plant_assimilation_over_distribution=TRUE)
  s1 <- Strategy(control=c1)
  p1 <- Plant(s1)

  c2 <- Control(plant_assimilation_adaptive=FALSE,
                plant_assimilation_over_distribution=TRUE)
  s2 <- Strategy(control=c2)
  p2 <- Plant(s2)

  p1$height <- 10.0
  p2$height <- p1$height
  env <- test_environment(p1$height)

  p1$compute_vars_phys(env)
  p2$compute_vars_phys(env)
  p1_phys <- p1$vars_phys
  p2_phys <- p2$vars_phys

  ## Result is similar but not identical:
  expect_that(p2_phys[["assimilation"]],
              equals(p1_phys[["assimilation"]], tolerance=1e-3))
  expect_that(p2_phys, not(is_identical_to(p1_phys)))
})

test_that("Ode interface", {
  p <- Plant(Strategy())
  expect_that(p$ode_size, equals(5))
  expect_that(p$ode_state,
              equals(c(p$height, p$mortality, p$fecundity,
                       p$heartwood_area, p$heartwood_mass)))

  env <- test_environment(p$height * 10)
  p$compute_vars_phys(env)
  expect_that(p$ode_state,
              equals(c(p$height, p$mortality, p$fecundity,
                       p$heartwood_area, p$heartwood_mass)))
  phys <- as.list(p$vars_phys)
  growth <- as.list(p$vars_growth)
  expect_that(p$ode_rates,
              equals(c(phys$height_growth_rate,
                       phys$mortality_rate,
                       phys$fecundity_rate,
                       growth$dheartwood_area_dt,
                       growth$dheartwood_mass_dt)))

  state_new <- c(p$height * 2, runif(p$ode_size - 1L))
  p$ode_state <- state_new
  expect_that(p$ode_state, is_identical_to(state_new))

  expect_that(p$ode_names,
              is_identical_to(c("height", "mortality", "fecundity",
                                "heartwood_area", "heartwood_mass")))
})
