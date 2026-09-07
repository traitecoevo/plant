# Guards the numbers TF24_flexible_allometry_demo.qmd asserts in prose.
#
# The demo's chunks all run, but a chunk that runs is not a chunk that is
# checked: prose beside a number survives the number being falsified. One claim
# in the first draft did exactly that (a Huber-value rise quoted from a shorter
# drought than the one plotted), which is why this file exists.

demo_helpers <- test_path("..", "..", "overstorey_staging",
                          "allometry_demo_helpers.R")

test_that("the flexible-allometry demo's claims still hold", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  history <- demo_soil_history(wet = 0.30, dry = 0.11,
                               years_wet = 2, years_dry = 4)
  fixed <- allometry_trajectory(0.0, history, dt = 0.25)
  flex  <- allometry_trajectory(1.0, history, dt = 0.25)

  ## Tolerances are loose because these are integrated trajectories, not
  ## closed forms; what is being guarded is the claim, not the digits.

  ## "thins to about 19 per cent of what its height prefers"
  expect_equal(min(flex$canopy_fraction), 0.19, tolerance = 0.15)
  ## the fixed plant does not move off its allometry at all
  expect_equal(max(abs(fixed$canopy_fraction - 1)), 0, tolerance = 1e-8)

  ## "rises to about 2.4 times its preferred value"
  expect_equal(max(flex$huber_ratio), 2.4, tolerance = 0.15)
  ## and psi >= 0 is invariant, so the ratio never drops below 1
  expect_gte(min(flex$huber_ratio), 1 - 1e-8)

  ## "rebuilds to ~98 per cent of preferred"
  expect_gt(tail(flex$canopy_fraction, 1), 0.9)

  ## Height never falls, in either model.
  expect_true(all(diff(flex$height) >= -1e-10))
  expect_true(all(diff(fixed$height) >= -1e-10))

  ## THE ISSUE THE DEMO EXISTS TO FLAG: thinning buys ~1 per cent survival and
  ## costs ~2 m of height. If either of these moves materially, the demo's
  ## conclusion has changed and its prose needs rewriting -- which is the point
  ## of asserting them.
  expect_equal(max(flex$survivorship / fixed$survivorship), 1.01,
               tolerance = 0.05)
  expect_equal(tail(fixed$height, 1) - tail(flex$height, 1), 2.07,
               tolerance = 0.4)
})

test_that("the demo's gate map brackets the replacement gate", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  gate <- allometry_gate_map(a_pl0 = 1.0, height = 5, n = 21L)

  ## Opposite signs: leaf area thins while sapwood per leaf area rises.
  expect_lt(gate$dphi[[1]], 0)
  expect_gt(gate$dpsi[[1]], 0)

  ## "by r = 0.33 the rates are already ~1e-7" -- the gate is sharp, which is
  ## the demo's explanation for why the response is small in most stands.
  near_third <- gate[which.min(abs(gate$r - 1 / 3)), ]
  expect_lt(abs(near_third$dphi), 1e-5)

  ## Ample reserves: no thinning at all.
  expect_equal(gate$dphi[[nrow(gate)]], 0, tolerance = 1e-12)

  ## Height growth is non-negative throughout.
  expect_true(all(gate$dheight >= 0))
})
