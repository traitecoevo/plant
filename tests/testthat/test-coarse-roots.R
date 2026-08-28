## Coarse (structural) root pool -- issue #349.
##
## FF16 and TF24 gained a second root pool, m_cr = a_cr1 * m_s: coarse roots as
## the below-ground continuation of the sapwood cylinder, so they inherit its
## height scaling and the structural root:shoot ratio is size-invariant. Fine
## roots (a_r1 * A_l) are unchanged.
##
## The pool ships OFF (a_cr1 = 0), which splits these tests in two. Half of them
## check it really is inert at the defaults -- the reason FF16 keeps agreeing
## with the published 2012 model in test-strategy-ff16-reference-comparison.R
## and no scientific version moved. The other half switch it on, because a
## default-off feature only ever run at its default is not tested at all.

eta_c_of <- function(p) 1 - 2 / (1 + p$eta) + 1 / (1 + 2 * p$eta)

## Published FF16 allometric derivatives, written out here rather than read back
## from the strategy, so the C++ is compared against the documented formulas
## (inst/docs/FF16/FF16-eqns.csv) and not against itself.
dmass_sapwood_darea_leaf_r <- function(p, area_leaf) {
  p$rho * eta_c_of(p) * p$a_l1 * p$theta * (p$a_l2 + 1) * area_leaf^p$a_l2
}

area_leaf_r <- function(p, height) (height / p$a_l1)^(1 / p$a_l2)

## Sum of d(component mass)/d(area_leaf) over the live pools -- the denominator
## of darea_leaf_dmass_live, i.e. the mass that must be built alongside each new
## unit of leaf area.
dmass_live_darea_leaf_r <- function(p, height) {
  dms <- dmass_sapwood_darea_leaf_r(p, area_leaf_r(p, height))
  p$lma + dms + p$a_b1 * dms + p$a_r1 + p$a_cr1 * dms
}

allometry_of <- function(type, s, height) {
  z <- rep(0, length(height))
  switch(type,
         FF16 = FF16_strategy_expand_allometry(s, height, z, z),
         TF24 = TF24_strategy_expand_allometry(s, height, z, z))
}

heights <- seq(0.5, 20, length.out = 40)

test_that("the coarse-root pool is absent at the shipped defaults", {
  for (type in c("FF16", "TF24")) {
    s <- get(paste0(type, "_Strategy"))()
    expect_identical(s$pars$a_cr1, 0, info = type)

    al <- allometry_of(type, s, heights)
    expect_true(all(al$mass_coarse_root == 0), info = type)
    ## Live and total mass therefore still read exactly as they did before the
    ## pool existed -- appending a term worth 0 changes no bits.
    expect_identical(al$mass_live,
                     al$mass_leaf + al$mass_sapwood + al$mass_bark + al$mass_root,
                     info = type)
  }
})

test_that("coarse-root mass is the sapwood fraction, and the mass budget closes", {
  a_cr1 <- 0.25
  for (type in c("FF16", "TF24")) {
    s <- get(paste0(type, "_Strategy"))()
    s$pars$a_cr1 <- a_cr1
    ## Non-zero heartwood, so mass_total and mass_above_ground are not just
    ## restatements of mass_live.
    mh <- seq(0, 5, length.out = length(heights))
    ah <- seq(0, 0.01, length.out = length(heights))
    al <- switch(type,
                 FF16 = FF16_strategy_expand_allometry(s, heights, ah, mh),
                 TF24 = TF24_strategy_expand_allometry(s, heights, ah, mh))

    expect_equal(al$mass_coarse_root, a_cr1 * al$mass_sapwood,
                 tolerance = 1e-14, info = type)
    ## Coarse roots are below ground: they belong to the total but not to the
    ## above-ground mass, and the two root pools are the whole difference.
    expect_equal(al$mass_total,
                 al$mass_above_ground + al$mass_root + al$mass_coarse_root,
                 tolerance = 1e-12, info = type)
    ## Size-invariant structural root:shoot -- the property that motivated
    ## scaling on sapwood rather than on leaf area.
    expect_equal(al$mass_coarse_root / al$mass_sapwood,
                 rep(a_cr1, length(heights)), tolerance = 1e-14, info = type)
  }
})

test_that("coarse roots slow height growth by exactly the allocation ratio", {
  ## Setting r_cr = k_cr = 0 makes the pool free to maintain, so net production
  ## is untouched and the ONLY route left is allocation: each new unit of leaf
  ## area now has to be built with a_cr1 * dm_s/dA_l of extra coarse root. The
  ## predicted slowdown is therefore the ratio of the two live-mass denominators,
  ## which is a closed form the strategy never computes.
  a_cr1 <- 0.3
  h0 <- 8

  dh_dt <- function(a_cr1, r_cr, k_cr) {
    s <- FF16_Strategy()
    s$pars$a_cr1 <- a_cr1
    s$pars$r_cr <- r_cr
    s$pars$k_cr <- k_cr
    env <- Environment("FF16")
    env$set_fixed_environment(1.0, height_max = 150)
    pl <- FF16_Individual(s)
    pl$set_state("height", h0)
    pl$compute_rates(env)
    list(dh = pl$rate("height"), P = pl$aux("net_mass_production_dt"), pars = s$pars)
  }

  off <- dh_dt(0,     0, 0)
  on  <- dh_dt(a_cr1, 0, 0)

  ## Free maintenance => identical carbon income.
  expect_identical(on$P, off$P)

  predicted <- dmass_live_darea_leaf_r(off$pars, h0) /
               dmass_live_darea_leaf_r(on$pars, h0)
  expect_equal(on$dh / off$dh, predicted, tolerance = 1e-12)
  expect_lt(on$dh, off$dh)

  ## The ratio alone would survive an error common to both denominators, so
  ## reconstruct dh/dt outright from the published chain
  ##     dH/dt = dH/dA_l * P * (1 - r(H)) / (dM_live/dA_l)
  ## and check the absolute value with the pool switched on.
  dheight_darea_leaf_r <- function(p, height) {
    p$a_l1 * p$a_l2 * area_leaf_r(p, height)^(p$a_l2 - 1)
  }
  frac_growth_r <- function(p, height) {
    1 - p$a_f1 / (1 + exp(p$a_f2 * (1 - height / p$hmat)))
  }
  expect_equal(on$dh,
               dheight_darea_leaf_r(on$pars, h0) * on$P *
                 frac_growth_r(on$pars, h0) /
                 dmass_live_darea_leaf_r(on$pars, h0),
               tolerance = 1e-12)
})

test_that("coarse roots also cost carbon to maintain", {
  ## At the shipped woody rates (r_cr = r_s, k_cr = k_s) the pool costs on both
  ## sides: respiration and turnover cut net production as well. So switching it
  ## on must slow growth strictly more than the allocation effect alone.
  h0 <- 8
  rates <- function(a_cr1, free) {
    s <- FF16_Strategy()
    s$pars$a_cr1 <- a_cr1
    if (free) { s$pars$r_cr <- 0; s$pars$k_cr <- 0 }
    env <- Environment("FF16")
    env$set_fixed_environment(1.0, height_max = 150)
    pl <- FF16_Individual(s)
    pl$set_state("height", h0)
    pl$compute_rates(env)
    c(dh = pl$rate("height"), P = pl$aux("net_mass_production_dt"))
  }

  off       <- rates(0.0, FALSE)
  alloc     <- rates(0.3, TRUE)
  alloc_mnt <- rates(0.3, FALSE)

  expect_lt(alloc_mnt[["P"]], off[["P"]])       # maintenance eats production
  expect_lt(alloc_mnt[["dh"]], alloc[["dh"]])   # ... on top of the allocation cost
  expect_lt(alloc[["dh"]], off[["dh"]])
})

test_that("coarse roots stay out of TF24's water uptake", {
  ## TF24 distributes FINE-root mass over soil depth to drive uptake, and reads
  ## a_r1 directly to do it. Coarse roots are structure, not absorbing surface,
  ## so at a fixed size every water-side quantity must be bit-identical whether
  ## the pool is on or off.
  water_aux <- function(a_cr1) {
    s <- TF24_Strategy()
    s$pars$a_cr1 <- a_cr1
    env <- Environment("TF24")
    env$set_fixed_environment(1.0, height_max = 150)
    env$set_soil_water_state(rep(0.4, env$get_soil_number_of_depths()))
    env$time <- 5
    ind <- Individual("TF24", "TF24_Env")(s)
    ind$set_state("height", 8)
    ind$compute_rates(env)
    vapply(c("root_mass", "transpiration", "E_up_", "opt_psi_stem",
             "opt_root_psi", "assimilation"),
           function(nm) ind$aux(nm), numeric(1))
  }

  expect_identical(water_aux(0.3), water_aux(0.0))
})
