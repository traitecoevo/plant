context("Refrence Comparison-FF16")

test_that("FF16_Strategy parameters agree with reference model", {
  cmp <- make_reference_plant("FF16")
  cmp_pars <- cmp$get_parameters()

  s <- FF16_Strategy()

  ## Expect that all parameters in the R version are found in the C++
  ## version, *except* for n_area
  v <- setdiff(names(cmp_pars), "n_area")
  expect_true(all(v %in% names(s)))

  ## And v.v., except for a few additions:
  extra <- c("control", "S_D")
  common <- setdiff(names(s), extra)
  expect_true(all(extra %in% names(s)))
  expect_true(all(common %in% names(cmp_pars)))

  ## The C++ version should have no NA values by this point.
  expect_false(any(sapply(s[common], is.na)))

  ## And neither should the R version.
  expect_false(any(sapply(cmp_pars, is.na)))


  ## And demand that all parameters agree.
  expect_equal(s[v], cmp_pars[v], tolerance=1e-13)
})

test_that("Reference comparison", {
  cmp <- make_reference_plant("FF16")
  s <- FF16_Strategy()
  p <- FF16_PlantPlus(s)

  expect_identical(p$strategy, s)

  ## Set the height to something (here 10)
  h0 <- 10
  p$height <- h0

  vars <- p$internals

  expect_identical(vars[["height"]], h0)
  expect_equal(vars[["area_leaf"]], cmp$LeafArea(h0))
  expect_equal(vars[["mass_leaf"]], cmp$LeafMass(cmp$traits$lma, cmp$LeafArea(h0)))
  expect_equal(vars[["mass_sapwood"]], cmp$SapwoodMass(cmp$traits$rho, cmp$LeafArea(h0), h0))
  expect_equal(vars[["mass_bark"]], cmp$BarkMass(cmp$traits$rho, cmp$LeafArea(h0), h0))
  expect_equal(vars[["mass_root"]], cmp$RootMass(cmp$LeafArea(h0)))
  expect_equal(vars[["mass_live"]], cmp$LiveMass(cmp$traits, cmp$LeafArea(h0)))
  expect_equal(vars[["area_bark"]], cmp$area_bark(h0))
  expect_equal(vars[["area_sapwood"]], cmp$area_sapwood(h0))
  expect_equal(vars[["area_stem"]], cmp$area_stem(h0))

  expect_identical(p$height, vars[["height"]])
  expect_identical(p$area_leaf, vars[["area_leaf"]])

  ## Heartwood function
  ## TODO: Check with Daniel here -- might need updating.
  ## Expect zero unless it has been set otherwise
  expect_identical(vars[["area_heartwood"]], 0.0)
  HA0 <- 1e-3
  p$area_heartwood <- HA0
  ## TODO: This is due to issues with how we think about size; see notes
  ## in plant.cpp around set_height and set_mass_heartwood.  Note that
  ## when running as an ODE, this gives the *wrong answer*.
  h <- p$height
  p$height <- h + .1 # trick plant into recomputing all size variables
  p$height <- h
  vars <- p$internals
  expect_identical(vars[["area_heartwood"]], HA0)
  expect_equal(vars[["area_stem"]], cmp$area_stem(h0) + HA0)
  # set heartwood back at zero for subsequent tests
  p$area_heartwood <- 0

  env <- test_environment(h0)
  light_env <- attr(env, "light_env") # underlying function

  ## The R model computes A_lf * area_leaf * a_y * a_bio, wheras we just
  ## compute A_lf; will have to correct some numbers.
  cmp_const <- s$a_y * s$a_bio

  ## TODO: Check what growth variables look like before running
  ## through with environment.  Most are NA_real_, but some are not

  ## Compute the physiological variables and retreive them.
  p$compute_vars_phys(env)
  p$compute_vars_growth() # NOTE: Compute immediately *after* vars_phys
  vars <- p$internals

  ## 1. Assimilation:
  cmp_assimilation_plant <- cmp$assimilation.plant(h0, light_env)
  expect_equal(vars[["assimilation"]], cmp_assimilation_plant / cmp_const)

  ## 2. Respiration:
  cmp_respiration <- cmp$respiration.given.height(cmp$traits, h0)
  expect_equal(vars[["respiration"]], cmp_respiration / cmp_const)

  ## 3. Turnover:
  cmp_turnover <- cmp$turnover.given.height(cmp$traits, h0)
  expect_equal(vars[["turnover"]], cmp_turnover)

  ## 4. Net production:
  cmp_net_mass_production_dt <- cmp$net.production(cmp$traits, h0, light_env)
  expect_equal(vars[["net_mass_production_dt"]], cmp_net_mass_production_dt, tolerance=1e-7)

  ## 5. Reproduction fraction
  cmp_fraction_allocation_reproduction <-
    cmp$ReproductiveAllocation(cmp$traits$hmat,h0)
  expect_equal(vars[["fraction_allocation_reproduction"]], cmp_fraction_allocation_reproduction)

  ## 6. Fecundity rate
  cmp_fecundity_dt <- cmp$fecundity_dt(cmp$traits, h0, light_env)
  expect_equal(vars[["fecundity_dt"]], cmp_fecundity_dt, tolerance=1e-7)

  ## 8. Growth rate for height
  cmp_height_dt <- cmp$height.growth.dt(cmp$traits, h0, light_env)
  expect_equal(vars[["height_dt"]], cmp_height_dt, tolerance=1e-7)

  cmp_height_dt <-
    cmp$height.growth.dt.via.area.leaf(cmp$traits, h0, light_env)
  expect_equal(vars[["height_dt"]], cmp_height_dt, tolerance=1e-7)

  ## 9. Mortality rate
  cmp_mortality_dt <- cmp$mortality.dt(cmp$traits, h0, light_env)
  expect_equal(vars[["mortality_dt"]], cmp_mortality_dt)

  ## 10. Archietcural layout
  cmp_dheight_darea_leaf <- cmp$dHdA(cmp$LeafArea(h0))
  expect_equal(vars[["dheight_darea_leaf"]], cmp_dheight_darea_leaf)

  ## 11. Sapwood mass per leaf mass
  cmp_dmass_sapwood_darea_leaf<- cmp$sapwood.per.leaf.area(cmp$traits, h0)
  expect_equal(vars[["dmass_sapwood_darea_leaf"]], cmp_dmass_sapwood_darea_leaf)

  ## 12. Bark mass per leaf mass
  cmp_dmass_bark_darea_leaf <- cmp$bark.per.leaf.area(cmp$traits, h0)
  expect_equal(vars[["dmass_bark_darea_leaf"]], cmp_dmass_bark_darea_leaf)

  ## 12. Root mass per leaf mass
  cmp_dmass_root_darea_leaf <- cmp$root.per.leaf.area(cmp$traits, h0)
  expect_equal(vars[["dmass_root_darea_leaf"]], cmp_dmass_root_darea_leaf)

  ## 13. Leaf area growth rate
  cmp_area_leaf_dt <- cmp$area_leaf_dt(cmp$traits, h0, light_env)
  expect_equal(vars[["area_leaf_dt"]], cmp_area_leaf_dt, tolerance=1e-7)

  ## 14. sapwood area growth rate
  cmp_area_sapwood_dt <- cmp$area_sapwood_dt(cmp$traits, h0, light_env)
  expect_equal(vars[["area_sapwood_dt"]], cmp_area_sapwood_dt, tolerance=1e-7)

  ## 15. bark area growth rate
  cmp_area_bark_dt <- cmp$area_bark_dt(cmp$traits, h0, light_env)
  expect_equal(vars[["area_bark_dt"]], cmp_area_bark_dt, tolerance=1e-7)

  ## 16. heartwood area growth rate
  cmp_area_heartwood_dt <- cmp$area_heartwood_dt(cmp$traits, h0, light_env)
  expect_equal(vars[["area_heartwood_dt"]], cmp_area_heartwood_dt, tolerance=1e-7)

  ## 17. basal area growth rate
  cmp_area_stem_dt <- cmp$area_stem_dt(cmp$traits, h0, light_env)
  expect_equal(vars[["area_stem_dt"]], cmp_area_stem_dt, tolerance=1e-7)

  ## 18. change in basal diam per basal area
  cmp_ddiameter_stem_darea_stem <- cmp$ddiameter_stem_darea_stem(cmp$area_stem(h0))
  expect_equal(vars[["ddiameter_stem_darea_stem"]], cmp_ddiameter_stem_darea_stem, tolerance=1e-7)

  ## 18. basal diam growth rate
  cmp_diameter_stem_dt <- cmp$diameter_stem_dt(cmp$traits, h0, light_env)
  expect_equal(vars[["diameter_stem_dt"]], cmp_diameter_stem_dt, tolerance=1e-7)

  ## Check that height decomposition multiplies out to give right
  ## answer
  cmp <- prod(unlist(vars[c("dheight_darea_leaf",
                            "darea_leaf_dmass_live",
                            "fraction_allocation_growth")]),
              vars[[c("net_mass_production_dt")]])
  expect_equal(vars[["height_dt"]], cmp, tolerance=1e-7)
})
