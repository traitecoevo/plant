# Check C++ implemnetation againsyt indenpendent implementation in R


test_that("FF16_Model parameters agree with reference model", {
  cmp <- make_reference_plant("FF16")
  cmp_pars <- cmp$get_parameters()

  s <- FF16_Model()
  ## Biological parameters now live in the nested `pars` sub-object.
  sp <- s$pars

  ## Expect that all parameters in the R version are found in the C++
  ## version, *except* for n_area
  v <- setdiff(names(cmp_pars), c("n_area"))
  expect_true(all(v %in% names(sp)))

  ## And v.v., except for a few additions present in pars but not the reference.
  ## The coarse-root trio (#349) is absent from the 2012 reference model by
  ## construction: at the shipped default a_cr1 = 0 the pool contributes nothing,
  ## so the two implementations still agree on every growth quantity below.
  extra <- c("S_D", "recruitment_decay", "a_cr1", "r_cr", "k_cr")
  common <- setdiff(names(sp), extra)
  expect_true(all(extra %in% names(sp)))
  expect_true(all(common %in% names(cmp_pars)))

  ## The C++ version should have no NA values by this point.
  expect_false(any(sapply(sp[common], is.na)))

  ## And neither should the R version.
  expect_false(any(sapply(cmp_pars, is.na)))


  ## And demand that all parameters agree.
  expect_equal(sp[v], cmp_pars[v], tolerance=1e-13)
})

test_that("Reference comparison", {
  cmp <- make_reference_plant("FF16")
  s <- FF16_Model()
  # skip('reference plant is plant plus')
  # p <- FF16_IndividualPlus(s)
  p <- FF16_Individual(s)

  expect_identical(p$strategy, s)

  ## Set the height to something (here 10)
  h0 <- 10
  p$set_state("height", h0)

  expect_identical(p$state("height"), h0)
  # testing set auxiliary state as well as area_leaf/competition_effect depends on height only
  expect_equal(p$aux("competition_effect"), cmp$LeafArea(h0))
  ## Sizes are not states; they are derived from height and heartwood by
  ## the model's own allometry.
  la <- cmp$LeafArea(h0)
  size <- FF16_strategy_expand_allometry(s, h0, 0, 0)
  expect_equal(size$mass_leaf, cmp$LeafMass(cmp$traits$lma, la))
  expect_equal(size$mass_sapwood, cmp$SapwoodMass(cmp$traits$rho, la, h0))
  expect_equal(size$mass_bark, cmp$BarkMass(cmp$traits$rho, la, h0))
  expect_equal(size$mass_root, cmp$RootMass(la))
  expect_equal(size$mass_live, cmp$LiveMass(cmp$traits, la))
  expect_equal(size$area_bark, cmp$area_bark(h0))
  expect_equal(size$area_sapwood, cmp$area_sapwood(h0))
  expect_equal(size$area_stem, cmp$area_stem(h0))

  ## Heartwood: a new plant has none, and a value set survives a change of
  ## height, which refreshes the height-dependent auxiliary values.
  expect_identical(p$state("area_heartwood"), 0.0)
  HA0 <- 1e-3
  p$set_state("area_heartwood", HA0)
  p$set_state("height", h0 + 0.1)
  p$set_state("height", h0)
  expect_identical(p$state("area_heartwood"), HA0)
  expect_equal(FF16_strategy_expand_allometry(s, h0, HA0, 0)$area_stem,
               cmp$area_stem(h0) + HA0)
  # set heartwood back at zero for subsequent tests
  p$set_state("area_heartwood", 0)

  env <- Environment("FF16")
  light_env <- function(x) {
    lapply(x, env$light_availability$get_value_at_height) |>
     as.numeric()
  }
  
  p$compute_rates(env)

  ## Assimilation, respiration, turnover and the allocation fractions are
  ## not exposed on their own; net production combines the first three, and
  ## the fecundity and height rates below carry the allocation.

  ## Net production:
  cmp_net_mass_production_dt <- cmp$net.production(cmp$traits, h0, light_env)
  expect_equal(p$aux("net_mass_production_dt"), cmp_net_mass_production_dt,
               tolerance=1e-7)

  ## Fecundity rate
  cmp_fecundity_dt <- cmp$fecundity_dt(cmp$traits, h0, light_env)
  expect_equal(p$rate("fecundity"), cmp_fecundity_dt, tolerance=1e-7)

  ## Growth rate for height
  cmp_height_dt <- cmp$height.growth.dt(cmp$traits, h0, light_env)
  expect_equal(p$rate("height"), cmp_height_dt, tolerance=1e-7)

  cmp_height_dt <-
    cmp$height.growth.dt.via.area.leaf(cmp$traits, h0, light_env)
  expect_equal(p$rate("height"), cmp_height_dt, tolerance=1e-7)

  ## Mortality rate
  cmp_mortality_dt <- cmp$mortality.dt(cmp$traits, h0, light_env)
  expect_equal(p$rate("mortality"), cmp_mortality_dt)

  ## Rates and derivatives of size: the model does not expose these, so
  ## take them from its own allometry by central differences in height and
  ## heartwood area, carried forward by the model's own rates.
  dh <- 1e-5
  dha <- 1e-8
  by_height <- function(v) {
    (FF16_strategy_expand_allometry(s, h0 + dh, 0, 0)[[v]] -
       FF16_strategy_expand_allometry(s, h0 - dh, 0, 0)[[v]]) / (2 * dh)
  }
  by_heartwood <- function(v) {
    (FF16_strategy_expand_allometry(s, h0, dha, 0)[[v]] -
       FF16_strategy_expand_allometry(s, h0, -dha, 0)[[v]]) / (2 * dha)
  }
  size_rate <- function(v) {
    by_height(v) * p$rate("height") +
      by_heartwood(v) * p$rate("area_heartwood")
  }

  ## Architectural layout and tissue per unit leaf area
  expect_equal(1 / by_height("area_leaf"), cmp$dHdA(la), tolerance=1e-6)
  expect_equal(by_height("mass_sapwood") / by_height("area_leaf"),
               cmp$sapwood.per.leaf.area(cmp$traits, h0), tolerance=1e-6)
  expect_equal(by_height("mass_bark") / by_height("area_leaf"),
               cmp$bark.per.leaf.area(cmp$traits, h0), tolerance=1e-6)
  expect_equal(by_height("mass_root") / by_height("area_leaf"),
               cmp$root.per.leaf.area(cmp$traits, h0), tolerance=1e-6)

  ## Area growth rates
  expect_equal(size_rate("area_leaf"),
               cmp$competition_effect_dt(cmp$traits, h0, light_env), tolerance=1e-6)
  expect_equal(size_rate("area_sapwood"),
               cmp$area_sapwood_dt(cmp$traits, h0, light_env), tolerance=1e-6)
  expect_equal(size_rate("area_bark"),
               cmp$area_bark_dt(cmp$traits, h0, light_env), tolerance=1e-6)
  cmp_area_heartwood_dt <- cmp$area_heartwood_dt(cmp$traits, h0, light_env)
  expect_equal(p$rate("area_heartwood"), cmp_area_heartwood_dt, tolerance=1e-7)
  expect_equal(size_rate("area_stem"),
               cmp$area_stem_dt(cmp$traits, h0, light_env), tolerance=1e-6)

  ## Stem diameter
  expect_equal(by_heartwood("diameter_stem") / by_heartwood("area_stem"),
               cmp$ddiameter_stem_darea_stem(cmp$area_stem(h0)), tolerance=1e-6)
  expect_equal(size_rate("diameter_stem"),
               cmp$diameter_stem_dt(cmp$traits, h0, light_env), tolerance=1e-6)

  ## Height growth decomposes as dH/dA * dA/dM * (1 - reproductive
  ## allocation) * net production.
  expect_equal(p$rate("height"),
               (1 - cmp$ReproductiveAllocation(cmp$traits$hmat, h0)) *
                 cmp_net_mass_production_dt / by_height("mass_live"),
               tolerance=1e-6)

  ## Whole-plant light compensation point
  cmp_lcp <- cmp$WPLCP(cmp$traits, h0)
  expect_equal(p$resource_compensation_point(), cmp_lcp, tolerance=1e-6)
})
