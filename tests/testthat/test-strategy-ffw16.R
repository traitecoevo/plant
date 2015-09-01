context("Strategy-FFW16")

test_that("Defaults", {
  expected <- list(
    a_l2     = 0.306,
    S_D   = 0.25,
    a_y      = 0.7,
    a_l1     = 5.44,
    a_r1     = 0.07,
    a_b1      = 0.17,
    r_b   = 8024 / 608,
    r_l   = 39.27 / 0.1978791,
    r_r   = 217,
    r_s   = 4012/608,
    a_f3  = 3.0*3.8e-5,
    a_bio  = 0.0245,
    d_I   = 0.01,
    a_dG1   = 5.5,
    a_dG2   = 20,
    a_p1   = 151.177775377968,
    a_p2   = 0.204716166503633,
    a_f1   = 1,
    a_f2   = 50,
    a_d0   = 0.1,
    eta    = 12,
    hmat   = 16.5958691,
    k_b    = 0.2,
    k_l   = 0.4565855,
    k_r    = 1,
    k_s   = 0.2,
    lma    = 0.1978791,
    rho    = 608,
    omega  = 3.8e-5,
    theta  = 4669,
    control = Control())

  keys <- sort(names(expected))

  s <- FFW16_Strategy()
  expect_that(s, is_a("FFW16_Strategy"))

  expect_that(sort(names(s)), is_identical_to(keys))
  expect_that(unclass(s)[keys], is_identical_to(expected[keys]))
})

test_that("FFW16_Strategy parameters agree with reference model", {
  cmp <- make_reference_plant("FFW16")
  cmp_pars <- cmp$get_parameters()

  s <- FFW16_Strategy()

  ## Expect that all parameters in the R version are found in the C++
  ## version, *except* for n_area
  v <- setdiff(names(cmp_pars), "n_area")
  expect_that(all(v %in% names(s)), is_true())

  ## And v.v., except for a few additions:
  extra <- c("control", "S_D")
  common <- setdiff(names(s), extra)
  expect_that(all(extra %in% names(s)), is_true())
  expect_that(all(common %in% names(cmp_pars)),
              is_true())

  ## The C++ version should have no NA values by this point.
  expect_that(any(sapply(s[common], is.na)),
              is_false())

  ## And neither should the R version.
  expect_that(any(sapply(cmp_pars, is.na)),
              is_false())


  ## And demand that all parameters agree.
  expect_that(s[v], equals(cmp_pars[v], tolerance=1e-13))
})

test_that("Reference comparison", {
  cmp <- make_reference_plant("FFW16")
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


test_that("FFW16_Strategy hyper-parameterisation", {

  s <- FFW16_Strategy()

  # lma
  lma <- c(0.1,1)
  ret <- FFW16_hyperpar(trait_matrix(lma, "lma"), s)
  expect_that(ret[, "lma"], equals(lma))
  expect_that(ret[, "k_l"], equals(c(1.46678,0.028600), tolerance=1e-5))
  expect_that(ret[, "r_l"], equals(c(392.70, 39.27), tolerance=1e-5))
  expect_that(colnames(ret), equals(c("lma","k_l","r_l")))

  # wood density
  rho <- c(200,300)
  ret <- FFW16_hyperpar(trait_matrix(rho, "rho"), s)
  expect_that(ret[, "rho"], equals(rho))
  expect_that(ret[, "r_s"], equals(c(20.06000,13.37333), tolerance=1e-5))
  expect_that(ret[, "r_b"], equals(2*ret[, "r_s"]))
  expect_that(colnames(ret), equals(c("rho","r_s","r_b")))

  # narea
  narea <- c(2E-3,2.3E-3)
  ret <- FFW16_hyperpar(trait_matrix(narea, "narea"), s)
  expect_that(ret[, "narea"], equals(narea))
  expect_that(ret[, "r_l"], equals(c(212.2508, 244.0884), tolerance=1e-5))
  expect_that(ret[, "a_p1"], equals(c(162.2592, 188.1549), tolerance=1e-5))
  expect_that(ret[, "a_p2"], equals(c(0.220904, 0.259173), tolerance=1e-5))
  expect_that(colnames(ret), equals(c("narea", "a_p1", "a_p2", "r_l")))

  # seed mass
  omega <- 3.8e-5*c(1,2,3)
  ret <- FFW16_hyperpar(trait_matrix(omega, "omega"), s)
  expect_that(ret[, "omega"], equals(omega))
  expect_that(ret[, "a_f3"], equals(3*omega))
  expect_that(colnames(ret), equals(c("omega","a_f3")))

  # Todo: need to example with empty trait matrix, just not sure how to trigger!

})

