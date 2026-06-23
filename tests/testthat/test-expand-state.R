context("expand_state allometry")

# Oracle: the historical hard-coded R allometry formulas that *_expand_state()
# used before they were replaced by calls into the strategy's own C++ functions
# (src/strategy_expand.cpp). Kept here so the C++ implementation stays locked to
# the documented formulas -- if the C++ ever drifts, this fails.
allometry_oracle <- function(s, height, area_heartwood, mass_heartwood) {
  p <- s$pars
  eta_c <- 1 - 2 / (1 + p$eta) + 1 / (1 + 2 * p$eta)
  area_leaf <- (height / p$a_l1)^(1.0 / p$a_l2)
  mass_leaf <- area_leaf * p$lma
  area_sapwood <- area_leaf * p$theta
  mass_sapwood <- area_sapwood * height * eta_c * p$rho
  area_bark <- p$a_b1 * area_leaf * p$theta
  mass_bark <- area_bark * height * eta_c * p$rho
  area_stem <- area_bark + area_sapwood + area_heartwood
  diameter_stem <- sqrt(4 * area_stem / pi)
  mass_root <- p$a_r1 * area_leaf
  mass_live <- mass_leaf + mass_sapwood + mass_bark + mass_root
  mass_total <- mass_leaf + mass_bark + mass_sapwood + mass_heartwood + mass_root
  mass_above_ground <- mass_leaf + mass_bark + mass_sapwood + mass_heartwood
  list(area_leaf = area_leaf, mass_leaf = mass_leaf, area_sapwood = area_sapwood,
       mass_sapwood = mass_sapwood, area_bark = area_bark, mass_bark = mass_bark,
       area_stem = area_stem, diameter_stem = diameter_stem, mass_root = mass_root,
       mass_live = mass_live, mass_total = mass_total,
       mass_above_ground = mass_above_ground)
}

check_allometry <- function(strategy_fn, cpp_fn) {
  s <- strategy_fn()
  height <- seq(0.5, 20, length.out = 50)
  # vary heartwood so area_stem / mass_total / mass_above_ground are exercised
  area_heartwood <- seq(0, 0.01, length.out = 50)
  mass_heartwood <- seq(0, 5, length.out = 50)

  cpp <- cpp_fn(s, height, area_heartwood, mass_heartwood)
  oracle <- allometry_oracle(s, height, area_heartwood, mass_heartwood)

  expect_setequal(names(cpp), names(oracle))
  for (nm in names(oracle)) {
    expect_equal(cpp[[nm]], oracle[[nm]], tolerance = 1e-12, info = nm)
  }
}

test_that("FF16 C++ allometry matches the historical R formulas", {
  check_allometry(FF16_Strategy, FF16_strategy_expand_allometry)
})

test_that("TF24 C++ allometry matches the historical R formulas", {
  check_allometry(TF24_Strategy, TF24_strategy_expand_allometry)
})

test_that("FF16_expand_state adds the expected derived columns", {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar, birth_rate = list(20))
  out <- run_scm(p, Environment("FF16"), Control(), collect = TRUE)
  es <- FF16_expand_state(out)

  derived <- c("area_leaf", "mass_leaf", "area_sapwood", "mass_sapwood",
               "area_bark", "mass_bark", "area_stem", "diameter_stem",
               "mass_root", "mass_live", "mass_total", "mass_above_ground")
  expect_true(all(derived %in% names(es$species)))

  # Re-derive with the oracle from the (already present) state columns and the
  # species' strategy, and confirm the pipeline output matches.
  s <- es$p$strategies[[1]]
  sp <- es$species
  oracle <- allometry_oracle(s, sp$height, sp$area_heartwood, sp$mass_heartwood)
  for (nm in names(oracle)) {
    expect_equal(sp[[nm]], oracle[[nm]], tolerance = 1e-12, info = nm)
  }
})
