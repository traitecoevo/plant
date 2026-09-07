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

  ## Above ~10 m the elasticity exceeds 1: leaf area is actively
  ## counterproductive, and removing leaves raises TOTAL assimilation.
  expect_lt(eta[[1]], 1)      # 5 m
  expect_gt(eta[[3]], 1)      # 15 m

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

  for (nm in c("wet", "moist")) {
    d <- ctl[[nm]]
    expect_equal(d$huber_ratio[which.min(abs(d$R_s))],
                 d$huber_ratio[which.max(d$dheight)])
  }

  ## The optimum moves with soil water, and is nowhere near the pipe model.
  peak <- vapply(ctl, function(d) d$huber_ratio[which.max(d$dheight)], numeric(1))
  expect_equal(unname(round(peak[["wet"]], 2)), 1.28)
  expect_equal(unname(round(peak[["moist"]], 2)), 1.35)
  expect_gt(peak[["dry"]], peak[["moist"]])       # drier wants MORE stem

  ## The wart the demo admits to: at the pipe-model ratio in dry soil the plant
  ## cannot grow, so the controller is switched off rather than pointing uphill.
  expect_lt(ctl$dry$dheight[[1]], 1e-8)
  expect_lt(abs(ctl$dry$R_s[[1]]), 1e-8)
  expect_gt(ctl$dry$R_s[[which.min(abs(ctl$dry$huber_ratio - 1.35))]], 0.05)
})
