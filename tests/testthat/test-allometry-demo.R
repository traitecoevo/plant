# Guards the numbers TF24_flexible_allometry_demo.qmd asserts in prose, and
# checks that its "The mathematics" appendix IS the model.
#
# The demo's chunks all run, but a chunk that runs is not a chunk that is
# checked: prose beside a number survives the number being falsified. Two claims
# in earlier drafts did exactly that -- a Huber-value rise quoted from a shorter
# drought than the one plotted, and a caption promising arrows that were never
# drawn -- and both rendered clean.

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

  ## Tolerances are loose because these are integrated trajectories, not closed
  ## forms; what is guarded is the claim, not the digits.

  ## The canopy thins substantially, and the fixed plant does not move at all.
  expect_lt(min(flex$canopy_fraction), 0.4)
  expect_equal(max(abs(fixed$canopy_fraction - 1)), 0, tolerance = 1e-8)

  ## Stem per leaf rises, and psi >= 0 is invariant so it never drops below 1.
  expect_gt(max(flex$huber_ratio), 1.5)
  expect_gte(min(flex$huber_ratio), 1 - 1e-8)

  ## The canopy floor holds: thinning stops at a_pl3 of the preferred canopy,
  ## which is what keeps leaf area in a range where the arithmetic means
  ## something.
  expect_gte(min(flex$canopy_fraction), TF24_Strategy()$pars$a_pl3 - 1e-8)

  ## It rebuilds afterwards.
  expect_gt(tail(flex$canopy_fraction, 1), 0.8)

  ## Height never falls, in either model.
  expect_true(all(diff(flex$height) >= -1e-10))
  expect_true(all(diff(fixed$height) >= -1e-10))

  ## THE ISSUE THE DEMO EXISTS TO FLAG: thinning buys almost no survival, and
  ## costs height. If either moves materially the demo's conclusion has changed
  ## and its prose needs rewriting -- which is the point of asserting them.
  expect_lt(max(flex$survivorship / fixed$survivorship), 1.10)
  expect_gt(tail(fixed$height, 1) - tail(flex$height, 1), 1.0)
})

test_that("the gate responds to soil water, not to the reserve pool", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  gate <- allometry_gate_map(a_pl0 = 1.0, height = 5,
                             theta = c(0.10, 0.13, 0.15, 0.20, 0.30))

  ## Reserves are held at a_st3 of capacity across the whole sweep, so a
  ## reserve-gated version could not respond anywhere on this axis.
  wet <- gate[gate$theta_soil >= 0.20, ]
  dry <- gate[gate$theta_soil <= 0.13, ]

  expect_true(all(wet$margin > 0))          # leaves pay for themselves
  expect_true(all(dry$margin < 0))          # they do not
  expect_true(all(abs(wet$dphi) < 1e-9))    # so no thinning when wet
  expect_true(all(dry$dphi < -0.1))         # and thinning when dry
  expect_true(all(dry$dpsi > 0))            # stem per leaf rises
  expect_true(all(gate$dheight >= 0))       # height never falls
})

test_that("suppressed plants in a wet stand thin; dominants do not", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  thin <- allometry_self_thinning(a_pl0 = 1.0, mpl = 15)
  d <- thin$cohorts[thin$cohorts$time == max(thin$cohorts$time), ]
  d <- d[order(-d$height), ]

  ## Shade is as much the point as drought, and nothing in the model knows about
  ## canopy position: the light gradient does all of it.
  expect_equal(d$canopy_fraction[[1]], 1, tolerance = 1e-6)   # the dominant
  expect_lt(min(d$canopy_fraction), 0.9)                      # someone thinned

  ## The floor holds here too -- this is the case where, without it, a
  ## suppressed cohort reached canopy 0.000 with stem-per-leaf 1.3e9.
  expect_gte(min(d$canopy_fraction), TF24_Strategy()$pars$a_pl3 - 1e-8)
  expect_true(all(is.finite(d$huber_ratio)))
  expect_lt(max(d$huber_ratio), 1e3)
})

# An independent R implementation of the equations written out in the demo's
# "The mathematics" appendix, checked against what the C++ actually computes.
#
# Written from the prose rather than from the source, so a divergence between
# the two shows up here. Constants not exposed to R are repeated with the member
# they mirror named: if one changes this fails, which is the intent.
tf24_departure_rates_r <- function(s, height, phi, psi, tracked, P, A) {
  p <- s$pars
  storage_gate_width    <- 0.1    # TF24_Strategy::storage_gate_width
  storage_prod_eps      <- 1e-4   # TF24_Strategy::storage_prod_eps
  plasticity_gate_width <- 0.02   # TF24_Strategy::plasticity_gate_width
  eta_c <- 1 - 2 / (1 + p$eta) + 1 / (1 + 2 * p$eta)   # CanopyShape::eta_c

  ## The gate reads the TRACKED marginal balance, and its effect tapers to
  ## nothing at the canopy floor.
  rho_carbon <- 1 - p$a_pl0 / (1 + exp((tracked - p$a_pl1) / plasticity_gate_width))
  shed_gate  <- 1 - exp((log(p$a_pl3) - phi) / p$a_pl2)
  rho_l <- 1 - (1 - rho_carbon) * shed_gate
  rho_s <- rho_l * exp(-psi / p$a_pl2)
  u_l   <- 1 - rho_l
  u_s   <- 1 - rho_s

  ## Reserves still gate growth, as they did before.
  A_s   <- p$theta * A * exp(psi)
  S_max <- p$a_st1 * A_s * height * eta_c * p$rho
  ## (r is supplied by the caller through P's own individual; recomputed here
  ## only for the growth gate.)
  list(rho_l = rho_l, u_l = u_l, u_s = u_s, S_max = S_max,
       finish = function(r) {
         G     <- 1 / (1 + exp(-(r - p$a_st2) / storage_gate_width))
         Ppos  <- 0.5 * (P + sqrt(P^2 + storage_prod_eps^2))
         F     <- Ppos * G
         f_g   <- 1 - p$a_f1 / (1 + exp(p$a_f2 * (1 - height / p$hmat)))
         sigma <- rho_l * (1 - exp(phi / p$a_pl2))
         c_reb <- p$lma + p$a_r1 +
           p$theta * height * eta_c * p$rho * (1 + p$a_b1)
         b <- sigma * F * f_g / (c_reb * A)
         list(dphi = b - u_l * p$k_l,
              dpsi = u_l * p$k_l - u_s * p$k_s + b * (exp(-psi) - 1))
       })
}

test_that("the demo's equations reproduce the C++ departure rates", {
  skip_on_cran()
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_fixed_environment(1.0, 40)

  s <- TF24_Strategy(collect_all_auxiliary = TRUE)
  s$pars$a_pl0 <- 1.0

  grid <- expand.grid(height = c(1, 5, 15),
                      phi = c(0, -0.4, -1.5),
                      psi = c(0, 0.3, 1.2),
                      theta_soil = c(0.11, 0.15, 0.30))

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    env$set_soil_water_state(rep(g$theta_soil, 5))
    ind <- TF24_Individual(s)
    ind$set_state("height", g$height)
    ind$set_state("log_area_leaf_departure", g$phi)
    ind$set_state("log_area_sapwood_departure", g$psi)
    ind$set_initial_states(env)
    ind$compute_rates(env)

    parts <- tf24_departure_rates_r(
      s, g$height, g$phi, g$psi,
      tracked = ind$state("leaf_marginal_return_tracked"),
      P = ind$aux("net_mass_production_dt"),
      A = ind$aux("competition_effect"))
    want <- parts$finish(ind$state("storage") / parts$S_max)

    label <- sprintf("h=%g phi=%g psi=%g soil=%g",
                     g$height, g$phi, g$psi, g$theta_soil)
    expect_equal(ind$rate("log_area_leaf_departure"), want$dphi,
                 tolerance = 1e-10, info = label)
    expect_equal(ind$rate("log_area_sapwood_departure"), want$dpsi,
                 tolerance = 1e-10, info = label)
  }
})

test_that("the shedding criterion's carbon decomposition is exact", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  ## P = c*abar*A - g(h)*A - s(h,A_s). This is the model's own budget
  ## rearranged, not an approximation of it, and the demo's criterion is
  ## nothing but its derivative -- so if this stops holding exactly, the whole
  ## analytical section is void rather than merely imprecise.
  for (h in c(5, 10, 15, 20)) {
    z <- allometry_carbon_terms(h, theta_soil = 0.20)
    expect_equal(z$gain - z$g - z$s, z$P, tolerance = 1e-10,
                 info = sprintf("height %g", h))
  }
})

test_that("the elasticity rises with height, and the criterion tracks dP/dA", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  z <- lapply(c(5, 10, 15, 20), allometry_carbon_terms, theta_soil = 0.20)
  eta   <- vapply(z, `[[`, numeric(1), "eta")
  kappa <- vapply(z, `[[`, numeric(1), "kappa")
  dPdA  <- vapply(z, `[[`, numeric(1), "dPdA")

  ## kmax ~ 1/h, so a taller plant is more supply-limited and relieving that
  ## limitation is worth more. This is what closes the gap in which tall plants
  ## died without ever shedding.
  expect_true(all(diff(eta) > 0))
  expect_true(all(diff(kappa) > 0))

  ## Somewhere above mid-height the elasticity passes 1, and beyond that leaf
  ## area is actively counterproductive: removing leaves raises TOTAL
  ## assimilation. Where the crossing sits is set by the height-resistance
  ## relation, so it moved when #617 replaced that with a path integral (eta at
  ## 15 m went from above 1 to 1.03, and it is 1.40 at 20 m). The crossing
  ## EXISTING is the claim; its location is the hydraulics' to set.
  expect_lt(eta[[1]], 1)                      # short plants: supply is not binding
  expect_gt(eta[[length(eta)]], 1)            # tall ones: it is

  ## The criterion and the measured derivative agree in sign. (Not independent
  ## -- eta is recovered through the decomposition, because the profit auxes
  ## cannot difference abar reliably -- but it does guard the arithmetic.)
  expect_equal(eta > 1 - kappa, dPdA < 0)

  ## A short 5 m plant in moist soil should not shed; a tall one should.
  expect_gt(dPdA[[1]], 0)
  expect_lt(dPdA[[4]], 0)
})

test_that("shedding is measured at fixed sapwood, which is what thinning does", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  ## Letting sapwood follow leaf area down the pipe model measures movement
  ## ALONG the allometry, not thinning, and gets the sign wrong where it
  ## matters. allometry_carbon_terms() holds A_s fixed and asserts it does; this
  ## checks the assertion is load-bearing by confirming the two differ.
  z <- allometry_carbon_terms(15, theta_soil = 0.20)
  expect_lt(z$dPdA, 0)                       # at fixed sapwood: shedding helps

  s <- TF24_Strategy(collect_all_auxiliary = TRUE)
  along <- function(phi) {
    e <- tf24_demo_env(0.20)
    i <- TF24_Individual(s)
    i$set_state("height", 15)
    i$set_state("log_area_leaf_departure", phi)   # psi left at 0: A_s follows A
    i$set_initial_states(e); i$compute_rates(e)
    c(A = i$aux("competition_effect"), P = i$aux("net_mass_production_dt"))
  }
  a <- along(0); b <- along(-0.02)
  expect_gt((b[["P"]] - a[["P"]]) / (b[["A"]] - a[["A"]]), 0)   # opposite sign
})

test_that("the demo's sapwood-controller claims hold", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  ## Guards the numbers the demo states in prose. A chunk that runs is not a
  ## chunk that is checked, and the claim here -- that the controller's zero
  ## sits on the growth peak -- is the one the whole derivation exists to make.
  ctl <- lapply(c(0.30, 0.20, 0.15), allometry_sapwood_controller)
  names(ctl) <- c("wet", "moist", "dry")

  ## ⚠️ THE ZERO IS FOUND BY SIGN CHANGE, not by min|R_s|. R_s has a SECOND zero
  ## at the thin-stem end, where supply collapses, production goes negative and
  ## the availability factor switches the controller off. That one is spurious
  ## for this purpose and it is also repelling (R_s jumps to +0.25 immediately
  ## above it), so a plant never rests there -- but min|R_s| finds it whenever
  ## the grid reaches far enough down, and reports the controller as broken.
  crossing <- function(d) {
    sg <- sign(d$R_s)
    i <- which(sg[-length(sg)] > 0 & sg[-1] <= 0)
    d$huber_ratio[i[length(i)]]
  }
  for (nm in names(ctl)) {
    d <- ctl[[nm]]
    expect_equal(crossing(d), d$huber_ratio[which.max(d$dheight)],
                 tolerance = 0.02)
  }

  ## The optimum moves with soil water: a drier plant wants MORE stem per leaf.
  ## Levels are not pinned -- they are the hydraulics' to set, and #617 moved
  ## them all -- but the ordering is the model's own claim.
  peak <- vapply(ctl, function(d) d$huber_ratio[which.max(d$dheight)], numeric(1))
  expect_lt(peak[["wet"]], peak[["moist"]])
  expect_lt(peak[["moist"]], peak[["dry"]])

  ## The wart the demo admits to: where the plant cannot grow at all, the
  ## controller is switched off rather than pointing uphill. After #617 a 10 m
  ## plant is solvent at soil 0.15, so this needs soil 0.12 to show -- which is
  ## itself the point, the band moved with the hydraulics.
  starved <- allometry_sapwood_controller(0.12)
  expect_true(all(starved$dheight < 1e-6))
  expect_true(all(abs(starved$R_s) < 1e-8))
})

test_that("the demo's 'does not buy survival' claims still hold", {
  skip_if_not(file.exists(demo_helpers), "demo helpers not present")
  skip_on_cran()
  source(demo_helpers, local = TRUE)

  ## This section rotted once already: it carried pre-#617 numbers describing a
  ## tall-plant cliff that no longer existed, and rendered clean the whole time
  ## because its tables are prose rather than chunks. These guard the claims.
  g <- allometry_drought_grid()

  ## 1. THE ORDERING, which is the section's decisive argument. At soil 0.145
  ##    every plant is healthy and none thins; by 0.140 a 10 m plant is at one
  ##    in a million and its canopy STILL has not moved. Mortality saturates at
  ##    a milder drought than shedding begins at.
  mild <- g[g$theta_dry == 0.145, ]
  expect_true(all(mild$min_canopy > 0.99))
  expect_true(all(abs(mild$gain - 1) < 0.01))

  mid10 <- g[g$theta_dry == 0.140 & g$height == 10, ]
  expect_equal(mid10$min_canopy, 1, tolerance = 1e-3)   # has not thinned
  expect_lt(mid10$S_fixed, 1e-5)                        # and is already dead

  ## 2. NO RESCUE. The large ratios sit between two very small numbers, so the
  ##    ratio alone must never be quoted as a survival benefit.
  expect_gt(max(g$gain), 5)                             # ratios do get large
  rescued <- g$S_fixed < 0.05 & g$S_shed > 0.2
  expect_false(any(rescued))

  ## 2b. The trajectory the section quotes: shedding ends BELOW parity, and the
  ##     height cost is about four metres. Measured on demo_soil_history's own
  ##     arguments -- quoting a number from a nearby-but-different run is how the
  ##     Huber-value claim went wrong in an earlier draft.
  hist <- demo_soil_history(wet = 0.30, dry = 0.11, years_wet = 2, years_dry = 4)
  tf <- allometry_trajectory(0.0, hist, dt = 0.25)
  tx <- allometry_trajectory(1.0, hist, dt = 0.25)
  ratio <- tail(tx$survivorship, 1) / tail(tf$survivorship, 1)
  expect_lt(ratio, 1)                                   # no gain; slightly worse
  expect_gt(ratio, 0.9)
  drop <- tail(tf$height, 1) - tail(tx$height, 1)
  expect_gt(drop, 3);  expect_lt(drop, 5)               # ~4 m shorter
  expect_lt(min(tx$canopy_fraction), 0.25)              # it really does thin

  ## 3. Thinning improves the per-leaf carbon balance but cannot close it --
  ##    the argument that REPLACED "thinning cuts income just as fast" when the
  ##    stem path integral landed. An interior optimum, and still negative.
  per_area <- vapply(c(0, -0.7, -1.6, -2.3, -3.0, -4.0), function(k) {
    s <- TF24_Strategy(collect_all_auxiliary = TRUE)
    ind <- TF24_Individual(s)
    ind$set_state("height", 10)
    ind$set_state("log_area_leaf_departure", k)
    ind$set_state("log_area_sapwood_departure", -k * 0.35)
    env <- tf24_demo_env(0.13)
    ind$set_initial_states(env)
    ind$compute_rates(env)
    ind$aux("net_mass_production_dt") / ind$aux("competition_effect")
  }, numeric(1))
  expect_true(all(per_area < 0))                        # never reaches balance
  expect_gt(per_area[[4]], per_area[[1]])               # but thinning helps
  expect_gt(which.max(per_area), 1)                     # interior optimum
  expect_lt(which.max(per_area), length(per_area))
})
