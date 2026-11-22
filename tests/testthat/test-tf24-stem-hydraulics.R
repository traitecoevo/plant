# Height-dependent stem hydraulic resistance.
#
# TF24 used to impose resistance strictly linear in height:
#   leaf_specific_conductance_max = K_s * theta / (height * eta_c)
# which over a seedling-to-mulga range (0.3 -> 8 m) is a ~29x resistance
# increase across exactly the size range where the demography happens. It was
# never derived -- it is what you get by assuming a uniform tube.
#
# The replacement integrates two within-plant anatomical profiles (conduit
# widening, and the Huber-value profile) along the flow path, so the height
# exponent is derived from measurable quantities rather than assumed. See
# plant/models/stem_hydraulics.h.
#
# Setting D_c, theta_c and L_tip all to zero collapses the path integral back to
# the linear model bit for bit. Several tests below exist to keep it that way.

# eta_c at the default eta = 12: the leaf-area-weighted mean height fraction of
# the Yokozawa crown, and therefore the fraction of full height that the
# representative flow path actually spans. Pinned here so that a change to the
# crown formula shows up as a stem-hydraulics failure too -- the two are coupled
# and nothing else asserts the coupling.
ETA_C_AT_12 <- 0.8861538461538462

test_that("stem path parameters carry their shipped defaults", {
  p <- TF24_Model()$pars
  # Widening is on; theta_c is zero, a uniform stem, so theta_tip is also the
  # whole-stem Huber value.
  expect_identical(p$D_c, 0.2)
  expect_identical(p$theta_c, 0)
  expect_identical(p$L_tip, 0.02)
  # K_s is now a TERMINAL-SEGMENT conductivity, back-derived from the old
  # whole-stem 1 so that resistance is unchanged at the anchor height, which the
  # path-integral test below checks.
  expect_equal(1 / p$K_s, 2.9781957045054700, tolerance = 1e-12)
})

test_that("zeroing the profile parameters recovers the height-linear model", {
  # Invariance criterion I7 survives the reparameterisation as an explicit
  # rather than as the default. expect_identical, not expect_equal: a value that
  # drifted to 1e-16 would still pass expect_equal while silently taking the
  # beta == 0 branch out of play, and the stand-level pins elsewhere sit at tolerance
  # 2e-2, which would absorb the consequence without anyone noticing.
  s <- TF24_Model()
  s$pars$D_c <- 0
  s$pars$theta_c <- 0
  s$pars$L_tip <- 0
  expect_identical(2 * s$pars$D_c + s$pars$theta_c, 0)
  expect_silent(TF24_Individual(s))

  eta_c <- 1 - 2 / (1 + s$pars$eta) + 1 / (1 + 2 * s$pars$eta)
  expect_identical(eta_c, ETA_C_AT_12)

  # 0.3941 is roughly TF24's birth height; 16.5958691 is its default hmat.
  for (h in c(0.2, 0.3941, 1, 5, 16.5958691, 40)) {
    L <- test_stem_effective_path_length(h * eta_c, 0, 0)
    # The path integral must return the OPERAND, not merely something equal to
    # it: no arithmetic at all is performed on the collapsed branch, which is
    # what removes any dependence on -ffp-contract fusing the caller's multiply.
    expect_identical(L, h * eta_c)
    # ...and the conductance the strategy forms from it must be, expression for
    # expression, the one the height-linear code computed.
    expect_identical(s$pars$K_s * s$pars$theta_tip / L,
                     s$pars$K_s * s$pars$theta_tip / (h * eta_c))
  }
})

test_that("a uniform stem builds the stem it always built, bit for bit", {
  # I7 for the carbon budget. At theta_c = 0 the Huber profile must be absent
  # from every structural quantity -- not merely close to absent -- in the
  # operation order the uniform stem always used, whatever L_tip is. The stand
  # pins elsewhere sit at tolerance 2e-2 and would not see a last-bit change.
  s <- TF24_Model()
  p <- s$pars
  expect_identical(p$theta_c, 0)
  # The strategy's own eta_c, so the comparison is of the same operands.
  eta_c <- canopy_profile(s, u = 1)$eta_c
  # Below and above the terminal segment (L_tip / eta_c = 0.0226 m).
  h <- c(0.01, 0.2, 0.3941, 1, 5, 16.5958691, 40)
  a <- plant:::TF24_strategy_expand_allometry(s, h, rep(0, length(h)),
                                              rep(0, length(h)))
  A <- (h / p$a_l1)^(1 / p$a_l2)
  expect_identical(a$area_sapwood, A * p$theta_tip)
  expect_identical(a$mass_sapwood, A * p$theta_tip * h * eta_c * p$rho)
  expect_identical(a$area_bark, p$a_b1 * A * p$theta_tip)
  expect_identical(a$mass_bark, p$a_b1 * A * p$theta_tip * h * eta_c * p$rho)
  expect_identical(a$huber_value_base, rep(p$theta_tip, length(h)))
  # And the hand-written marginal cost, at a resting and a displaced departure.
  for (psi in c(0, 0.3)) {
    for (A1 in A) {
      st <- plant:::test_tf24_stem_profile(s, A1, psi)
      expect_identical(st$dmass_sapwood_darea_leaf,
                       p$rho * eta_c * p$a_l1 * (p$theta_tip * exp(psi)) *
                         (p$a_l2 + 1) * A1^p$a_l2)
    }
  }
})

test_that("resistance is unchanged at the anchor height and rotates about it", {
  # The single point of agreement between the old and new models. R_L depends on
  # theta and K_s only through their ratio, so there is exactly one free scalar
  # and agreement holds at exactly ONE height -- never two. Below the anchor
  # every plant is more resistant than under the old model, above it every plant
  # is less: the reparameterisation is a rotation, not a rescaling.
  p <- TF24_Model()$pars
  eta_c <- canopy_profile(TF24_Model(), u = 1)$eta_c
  beta <- 2 * p$D_c + p$theta_c

  # theta cancels; K_s_old was 1.
  r_old <- function(h) test_stem_effective_path_length(h * eta_c, p$L_tip, 0) / 1
  r_new <- function(h) test_stem_effective_path_length(h * eta_c, p$L_tip, beta) / p$K_s

  expect_equal(r_new(TF24_H_ANCHOR) / r_old(TF24_H_ANCHOR), 1, tolerance = 1e-12)

  # The rotation, pinned. Anchored at 1 m, only sub-metre plants pay more than
  # they did; everything taller pays progressively less, down to 0.21x at 60 m.
  # That is the whole behavioural content of the change and it must not move
  # unnoticed.
  expect_equal(r_new(0.3941) / r_old(0.3941), 1.3757, tolerance = 1e-3)
  expect_equal(r_new(8) / r_old(8), 0.4615, tolerance = 1e-3)
  expect_equal(r_new(60) / r_old(60), 0.2100, tolerance = 1e-3)

  # Monotone through the anchor: strictly more resistant below, less above.
  expect_gt(r_new(0.5) / r_old(0.5), 1)
  expect_lt(r_new(30) / r_old(30), 1)
})

test_that("the K_s reparameterisation leaves the vulnerability curve alone", {
  # Invariance criterion I12. make_TF24_hyperpar derives the whole vulnerability
  # curve from K_s via stem_P50 = 10^(B_Hv1 + B_Hv2*log10(K_s)). Feeding a
  # terminal-segment K_s through the un-re-anchored relation would have moved
  # stem_P50 from 2.889 to 3.593 MPa -- buying a height exponent and silently
  # selling the safety margin of a model whose whole subject is hydraulic
  # limitation. B_Hv1 was shifted to hold it.
  s <- TF24_Model()
  m <- trait_matrix(s$pars$lma, "lma")
  derived <- default_hyperpar("TF24")(m, s, filter = FALSE)
  expect_equal(derived[, "stem_P50"][[1]], 2.888725665336019, tolerance = 1e-10)

  # And the un-re-anchored value, recorded so the size of the averted leak stays
  # visible: this is what stem_P50 would be at B_Hv1 = 0.4607063.
  expect_equal(10^(0.4607063 - 0.2 * log10(s$pars$K_s)), 3.5933287036419,
               tolerance = 1e-10)
})

test_that("stem path parameters are settable", {
  s <- TF24_Model()
  s$pars$D_c <- 0.3
  s$pars$L_tip <- 0.05
  expect_equal(s$pars$D_c, 0.3)
  expect_equal(s$pars$L_tip, 0.05)
})

test_that("the stem's areas and masses follow the Huber profile", {
  # theta(L) = theta_tip * (max(L, L_tip)/L_tip)^theta_c with L measured down
  # from the mean leaf height, the path the resistance integrates. Basal areas
  # are read at L = H*eta_c; masses are the profile's volume. Checked against
  # the profile integrated by R's quadrature, so the closed form is not checked
  # against itself, at exponents either side of the log limit theta_c = -1 and
  # on both sides of the terminal segment.
  base <- TF24_Model()
  p <- base$pars
  eta_c <- ETA_C_AT_12
  h <- c(0.01, 0.0225, 0.0227, 0.3941, 1, 16.5958691, 40)
  A <- (h / p$a_l1)^(1 / p$a_l2)
  for (theta_c in c(0.4, 1, 1.5, -0.3, -1)) {
    s <- base
    s$pars$theta_c <- theta_c
    a <- plant:::TF24_strategy_expand_allometry(s, h, rep(0, length(h)),
                                                rep(0, length(h)))
    L_top <- h * eta_c
    m <- ifelse(L_top > p$L_tip, (L_top / p$L_tip)^theta_c, 1)
    V <- vapply(L_top, function(Lt) {
      if (Lt <= p$L_tip) return(Lt)
      p$L_tip + stats::integrate(function(L) (L / p$L_tip)^theta_c, p$L_tip, Lt,
                                 rel.tol = 1e-12)$value
    }, numeric(1))
    info <- sprintf("theta_c = %g", theta_c)
    expect_equal(a$huber_value_base, p$theta_tip * m, tolerance = 1e-13, info = info)
    expect_equal(a$area_sapwood, A * p$theta_tip * m, tolerance = 1e-13, info = info)
    expect_equal(a$area_bark, p$a_b1 * A * p$theta_tip * m, tolerance = 1e-13,
                 info = info)
    expect_equal(a$mass_sapwood, p$rho * A * p$theta_tip * V, tolerance = 1e-10,
                 info = info)
    expect_equal(a$mass_bark, p$a_b1 * p$rho * A * p$theta_tip * V,
                 tolerance = 1e-10, info = info)
  }
})

test_that("a positive theta_c tapers the stem: more sapwood at the base", {
  # The direction of the compensation: basal Huber value rises with height,
  # and the sapwood a tall plant pays for rises with it.
  s <- TF24_Model()
  s$pars$theta_c <- 0.4
  h <- c(1, 5, 16.5958691, 40)
  zero <- rep(0, length(h))
  tapered <- plant:::TF24_strategy_expand_allometry(s, h, zero, zero)
  uniform <- plant:::TF24_strategy_expand_allometry(TF24_Model(), h, zero, zero)
  expect_true(all(diff(tapered$huber_value_base) > 0))
  expect_true(all(tapered$mass_sapwood > uniform$mass_sapwood))
  expect_true(all(diff(tapered$mass_sapwood / uniform$mass_sapwood) > 0))
})

test_that("the hand-written sapwood cost is the derivative of the sapwood mass", {
  # dmass_sapwood_darea_leaf is written out rather than differentiated, so
  # nothing else fails if it drifts from mass_sapwood. Checked against a tangent
  # derivative of the mass along the resting allometry, H = a_l1 * A^a_l2, at
  # displaced departures and across the terminal segment.
  base <- TF24_Model()
  p <- base$pars
  h <- c(0.01, 0.3941, 1, 16.5958691, 40)
  A <- (h / p$a_l1)^(1 / p$a_l2)
  for (theta_c in c(0, 0.4, 1, 1.5, -0.3)) {
    s <- base
    s$pars$theta_c <- theta_c
    for (psi in c(0, 0.3)) {
      for (A1 in A) {
        st <- plant:::test_tf24_stem_profile(s, A1, psi)
        expect_equal(st$dmass_sapwood_darea_leaf,
                     st$dmass_sapwood_darea_leaf_tangent, tolerance = 1e-12,
                     info = sprintf("theta_c = %g, psi = %g, A = %g",
                                    theta_c, psi, A1))
      }
    }
  }
})

test_that("theta_c's derivative is carried at a uniform stem and off it", {
  # At theta_c = 0 the value is the uniform stem's expression and the theta_c
  # column is attached by hand, as it is at theta_c = -1, the volume's log limit;
  # elsewhere the general closed forms carry it. All checked against a central
  # difference of the double values in theta_c.
  base <- TF24_Model()
  p <- base$pars
  h <- c(0.3941, 5, 16.5958691)
  A <- (h / p$a_l1)^(1 / p$a_l2)
  fields <- c("area_sapwood", "mass_sapwood", "area_bark", "mass_bark",
              "mass_live_per_area_leaf_at_height", "dmass_sapwood_darea_leaf")
  at <- function(theta_c, A1, psi) {
    s <- base
    s$pars$theta_c <- theta_c
    plant:::test_tf24_stem_profile(s, A1, psi)
  }
  d <- 1e-6
  for (theta_c in c(0, 0.4, -1)) {
    for (A1 in A) {
      for (psi in c(0, 0.3)) {
        ad <- at(theta_c, A1, psi)$by_theta_c
        up <- at(theta_c + d, A1, psi)
        dn <- at(theta_c - d, A1, psi)
        for (f in fields) {
          fd <- (up[[f]] - dn[[f]]) / (2 * d)
          expect_equal(ad[[f]], fd, tolerance = 1e-7,
                       info = sprintf("%s, theta_c = %g, A = %g, psi = %g",
                                      f, theta_c, A1, psi))
        }
      }
    }
  }
})

test_that("theta_c's stand-gradient column is the pinned-step difference", {
  # At the default theta_c = 0, where the column is attached by hand. The
  # difference itself scatters by about 1e-5 relative across steps of 5e-6 to
  # 2e-5 on this stand (net reproduction ratio -93.8150, -93.8141, -93.8151),
  # and the gradient sits inside that scatter, so 3e-5 is the difference's
  # floor rather than a loosened tolerance. A missing or wrong attached term
  # is out by order one.
  p <- tf24_fitness_stand_parameters()
  solver <- tf24_fitness_stand_run(p)
  metrics <- c("leaf_area", "1.net_reproduction_ratio")
  g <- stand_gradient(solver, metrics = metrics, traits = "1.theta_c",
                      hyperpar = NULL)
  expect_false(any(stand_gradient_refused(g)))
  pinned <- pin_steps(solver)
  measure <- function(q) {
    sv <- tf24_fitness_stand_run(q, record = FALSE)
    c(stand_census(sv), stand_fitness(sv))
  }
  at <- function(v) {
    q <- pinned
    q$strategies[[1]]$pars$theta_c <- v
    measure(q)
  }
  h <- 1e-5
  up <- at(h)
  dn <- at(-h)
  for (m in metrics) {
    fd <- (up[[m]] - dn[[m]]) / (2 * h)
    expect_lt(abs(g$gradient[m, "1.theta_c"] - fd), 3e-5 * abs(fd), label = m)
  }
})

test_that("a profiled stem builds, runs and needs a terminal segment", {
  s <- TF24_Model()
  s$pars$theta_c <- 0.4
  expect_silent(TF24_Individual(s))
  expect_true(is.finite(tf24_probe(s, rep(0.25, 5), height = 5)[["assimilation"]]))

  # theta_c = -2*D_c makes beta zero, so the path-length check does not ask for
  # a terminal segment, but the stem's profile still needs one.
  s2 <- TF24_Model()
  s2$pars$theta_c <- -0.4
  s2$pars$L_tip <- 0
  expect_identical(2 * s2$pars$D_c + s2$pars$theta_c, 0)
  expect_error(TF24_Individual(s2), "theta_c != 0 needs a terminal segment")

  s3 <- TF24_Model()
  s3$pars$theta_c <- NaN
  expect_error(TF24_Individual(s3), "theta_c must be finite")
})

test_that("theta was renamed theta_tip, and the old name is refused", {
  # Invariance criterion I11: a measured whole-stem or basal Huber value handed
  # over under the old name must not be read as a tip value.
  s <- TF24_Model()
  expect_null(s$pars$theta)

  # Beside theta_tip: old code setting a field nothing reads.
  s$pars$theta <- 1e-3
  expect_error(TF24_Individual(s), "`theta` is now `theta_tip`")

  # Alone, with theta_c = 0: an object saved before the rename, translated.
  old <- TF24_Model()
  pars <- old$pars
  names(pars)[names(pars) == "theta_tip"] <- "theta"
  pars$theta <- 3e-4
  old$pars <- pars
  expect_identical(TF24_Individual(old)$strategy$pars$theta_tip, 3e-4)

  # Alone, with theta_c != 0: tip and whole stem differ, so it is refused.
  pars$theta_c <- 0.4
  old$pars <- pars
  expect_error(TF24_Individual(old), "give theta_tip explicitly")

  # As a trait, named for what it became rather than its nearest spelling.
  p <- testing_parameters("TF24")
  expect_error(generate_strategy(p, trait_matrix(3e-4, "theta"), hyperpar = NULL),
               "theta \\(renamed theta_tip\\)")
})

test_that("the closed form matches numerical quadrature of the integrand", {
  # L_eff = INT_{L_tip}^{L_top} (L/L_tip)^(-beta) dL, checked against R's own
  # adaptive quadrature. beta > 1 is in the grid deliberately: that is the
  # saturating regime, where the exponent flips sign and both numerator and
  # denominator of the expm1 form go negative. It has no branch of its own, so
  # nothing else would catch a sign error there.
  quad <- function(L_top, L_tip, beta) {
    stats::integrate(function(L) (L / L_tip)^(-beta), L_tip, L_top,
                     rel.tol = .Machine$double.eps^0.75)$value
  }
  grid <- expand.grid(L_top = c(0.2, 1, 7, 40, 400),
                      L_tip = c(0.002, 0.02, 0.2),
                      beta = c(0.05, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0))
  grid <- grid[grid$L_tip < grid$L_top, ]
  for (i in seq_len(nrow(grid))) {
    with(grid[i, ],
         expect_equal(test_stem_effective_path_length(L_top, L_tip, beta),
                      quad(L_top, L_tip, beta),
                      tolerance = 1e-10,
                      info = sprintf("L_top=%g L_tip=%g beta=%g",
                                     L_top, L_tip, beta)))
  }
})

test_that("beta = 1 takes the logarithmic limit, continuously", {
  L_top <- 7
  L_tip <- 0.02
  f <- function(beta) test_stem_effective_path_length(L_top, L_tip, beta)

  expect_equal(f(1), L_tip * log(L_top / L_tip), tolerance = 1e-14)

  # The guarded branch must not be a discontinuity. This is the numerical face
  # of invariance criterion I10 (R_L continuously differentiable in H): a kink
  # in the parameter direction would generate singular strategies at the kink
  # just as mechanically as one in the height direction.
  for (d in c(1e-6, 1e-9, 1e-12)) {
    expect_equal(f(1 - d), f(1), tolerance = 10 * d)
    expect_equal(f(1 + d), f(1), tolerance = 10 * d)
  }

  # And this is why the implementation departs from the algebraic form in sec.
  # 4.4 of the design note. Written as L_tip^beta*(L_top^(1-beta) -
  # L_tip^(1-beta))/(1-beta), both powers approach 1 as beta -> 1 and
  # differencing them keeps only ~4 significant digits at 1-beta = 1e-12; the
  # expm1 form keeps all of them. If someone "simplifies" back to the note's
  # version, this assertion is what fails.
  expect_equal(f(1 - 1e-12), f(1), tolerance = 1e-11)
})

test_that("the path length's beta derivative is carried at the limits too", {
  # beta == 0 and beta == 1 take arms that return values not reading beta, so an
  # active scalar attaches the column by hand there. Without it the D_c and
  # theta_c columns lose the path's whole share at exactly those exponents and
  # read as answers. Checked against a central difference, which crosses into
  # the expm1 form on both sides.
  d <- 1e-6
  for (L_top in c(0.5, 7, 40)) {
    for (beta in c(0, 1, 0.4, 1.5)) {
      fd <- (test_stem_effective_path_length(L_top, 0.02, beta + d) -
               test_stem_effective_path_length(L_top, 0.02, beta - d)) / (2 * d)
      expect_equal(plant:::test_stem_effective_path_length_dbeta(L_top, 0.02, beta),
                   fd, tolerance = 1e-7,
                   info = sprintf("L_top = %g, beta = %g", L_top, beta))
    }
  }
})

test_that("path length is monotone and its height elasticity is the analytic one", {
  # The exact elasticity of the closed form, with x = L_top/L_tip and e = 1-beta:
  #   d ln L_eff / d ln L_top  =  e / (1 - x^-e)      (e != 0)
  #                            =  1 / log(x)          (e == 0, the log limit)
  #
  # Note what this says about sec. 4.3. The asymptotic exponent 1-beta is
  # approached only from ABOVE and only when e > 0; at x = 5000 and beta = 0.4
  # the elasticity is still 0.604, not 0.600. And for beta >= 1 it does not
  # approach 1-beta at all -- it tends to ZERO, because resistance approaches a
  # finite asymptote no matter how tall the plant grows. That is the saturating
  # regime, and an asymptotic reading of the table would have it backwards.
  elasticity_exact <- function(L_top, L_tip, beta) {
    x <- L_top / L_tip
    e <- 1 - beta
    if (e == 0) 1 / log(x) else e / (1 - x^(-e))
  }

  L_tip <- 0.02
  H <- 10^seq(-0.5, 2, length.out = 40)
  midpoints <- sqrt(H[-1] * H[-length(H)])   # geometric, matching diff(log(H))

  for (beta in c(0, 0.2, 0.4, 0.8, 1.0, 1.5)) {
    L <- vapply(H, function(h) test_stem_effective_path_length(h, L_tip, beta),
                numeric(1))
    expect_true(all(diff(L) > 0),
                info = sprintf("not increasing in H at beta=%g", beta))

    elasticity <- diff(log(L)) / diff(log(H))
    expect_true(all(diff(elasticity) <= 1e-9),
                info = sprintf("elasticity not non-increasing at beta=%g", beta))

    expected <- vapply(midpoints, elasticity_exact, numeric(1),
                       L_tip = L_tip, beta = beta)
    expect_equal(elasticity, expected, tolerance = 1e-3,
                 info = sprintf("elasticity off the analytic curve at beta=%g",
                                beta))
  }

  # The two limiting regimes, stated separately so a regression says which one
  # broke. Below beta = 1 the elasticity settles just above 1-beta; at and above
  # it, resistance saturates and the elasticity collapses towards zero.
  #
  # The approach is slow, and gets slower as beta grows: at an absurd H = 1e6 m
  # the elasticity is 2e-8 relative above its asymptote at beta = 0, but still
  # 3.0% above it at beta = 0.8. That is worth knowing before quoting the
  # asymptotic exponent for the recommended configuration -- over a real height
  # range the realised exponent is meaningfully steeper than 0.2.
  tall <- 1e6
  for (beta in c(0, 0.2, 0.4, 0.8)) {
    e_tall <- elasticity_exact(tall, L_tip, beta)
    expect_gt(e_tall, 1 - beta)
    expect_lt(e_tall / (1 - beta) - 1, 0.05)
  }
  for (beta in c(1.0, 1.5, 2.0)) {
    expect_lt(elasticity_exact(tall, L_tip, beta), 0.06)
  }

  # ...and decreasing in beta at fixed height: more compensation, less path.
  at_10m <- vapply(seq(0, 1.5, by = 0.05),
                   function(b) test_stem_effective_path_length(10, L_tip, b),
                   numeric(1))
  expect_true(all(diff(at_10m) < 0))
})

test_that("the conduit widening the law implies is reported at the ground", {
  # The diagnostic of #615's Phase 4: D_base / D_tip = (H/L_tip)^D_c, the
  # widening a plant of height H needs at the ground, L = H from the apex, to
  # set against measured vessels. It is 1 for a plant no taller than its
  # terminal segment and with widening off. Checked over L_tip, D_c and the
  # crown shape, so a hard-coded segment length, or the mean-path base
  # H*eta_c the resistance uses, fails here.
  h <- c(0.01, 0.0201, 0.3941, 1, 16.5958691, 60)
  zero <- rep(0, length(h))
  grid <- expand.grid(D_c = c(0.2, 0.3, 0), L_tip = c(0.02, 0.05), eta = c(12, 1))
  for (k in seq_len(nrow(grid))) {
    g <- grid[k, ]
    s <- TF24_Model()
    s$pars$D_c <- g$D_c
    s$pars$L_tip <- g$L_tip
    s$pars$eta <- g$eta
    a <- plant:::TF24_strategy_expand_allometry(s, h, zero, zero)
    expected <- ifelse(h > g$L_tip, (h / g$L_tip)^g$D_c, 1)
    expect_equal(a$conduit_widening_base, expected, tolerance = 1e-14,
                 info = sprintf("D_c = %g, L_tip = %g, eta = %g", g$D_c, g$L_tip, g$eta))
  }
  # Exactly 1 with widening off, and for a plant within its terminal segment.
  s <- TF24_Model()
  s$pars$D_c <- 0
  expect_identical(
    plant:::TF24_strategy_expand_allometry(s, h, zero, zero)$conduit_widening_base,
    rep(1, length(h)))
  expect_identical(
    plant:::TF24_strategy_expand_allometry(TF24_Model(), 0.01, 0, 0)$conduit_widening_base,
    1)
  # A missing height is not a plant within its terminal segment.
  expect_true(is.nan(
    plant:::TF24_strategy_expand_allometry(TF24_Model(), NaN, 0, 0)$conduit_widening_base))

  # At the defaults: the widening at 1 m, at hmat and at 60 m. A change to D_c
  # or L_tip that moves these moves the comparison with measured vessels.
  a <- plant:::TF24_strategy_expand_allometry(TF24_Model(), c(1, 16.5958691, 60),
                                              c(0, 0, 0), c(0, 0, 0))
  expect_equal(a$conduit_widening_base, c(2.1867, 3.8353, 4.9593), tolerance = 1e-4)

  # The subclasses report it too.
  a <- plant:::TF24f_strategy_expand_allometry(TF24f_Model(), 16.5958691, 0, 0)
  expect_equal(a$conduit_widening_base, (16.5958691 / 0.02)^0.2, tolerance = 1e-14)

  # It reaches a run's tidied output per cohort, through TF24_expand_state().
  p <- tf24_fitness_stand_parameters()
  res <- run(p, tf24_half_wet_environment(),
             Control(node_density_coordinate = "birth_date", ode_tol_rel = 1e-6,
                     ode_tol_abs = 1e-6), record = "steps")
  sp <- TF24_expand_state(tidy(res))$species
  ok <- is.finite(sp$height)
  expect_true(any(ok))
  expect_equal(sp$conduit_widening_base[ok],
               ifelse(sp$height[ok] > 0.02, (sp$height[ok] / 0.02)^0.2, 1),
               tolerance = 1e-12)
})

test_that("theta_c enters with the sign that compensates", {
  # theta(L) = theta_tip * (L/L_tip)^theta_c: the Huber value, sapwood area per
  # leaf area, RISES basipetally, so resistance per unit length falls as
  # (L/L_tip)^-(2*D_c + theta_c). A positive theta_c must therefore reduce the
  # effective path length, i.e. reduce resistance. Getting this sign backwards
  # inverts the compensation the whole exercise is about, and every other test
  # in this file would still pass -- they are all symmetric in beta.
  L_top <- 10
  L_tip <- 0.02
  none <- test_stem_effective_path_length(L_top, L_tip, 0)
  huber_only <- test_stem_effective_path_length(L_top, L_tip, 0.4)   # theta_c=0.4
  widen_only <- test_stem_effective_path_length(L_top, L_tip, 2 * 0.2)  # D_c=0.2
  expect_lt(huber_only, none)
  # ...and the two mechanisms are interchangeable at equal beta: 2*D_c and
  # theta_c enter the integral only through their sum.
  expect_identical(huber_only, widen_only)
})

test_that("a non-zero profile requires a terminal segment", {
  # The default now HAS widening on, so drop L_tip rather than raise D_c.
  s <- TF24_Model()
  s$pars$L_tip <- 0
  expect_error(TF24_Individual(s), "L_tip")

  s$pars$L_tip <- 0.02
  expect_silent(TF24_Individual(s))

  # Longer than the birth-size flow path: the plant would be shorter than one
  # terminal segment, the path length would go negative, and the symptom
  # without this guard is an unattributable NaN inside the leaf solver.
  s$pars$L_tip <- 100
  expect_error(TF24_Individual(s), "terminal segment")

  # The segment's length is checked at any profile exponent: at beta == 0 the
  # path is L_top - L_tip, so a long segment makes it negative there too.
  linear <- TF24_Model()
  linear$pars$D_c <- 0
  linear$pars$L_tip <- 100
  expect_error(TF24_Individual(linear), "terminal segment")
  linear$pars$L_tip <- 0.02
  expect_silent(TF24_Individual(linear))
})

test_that("a plant below its terminal segment is a domain error at run time", {
  # The seed is checked once, above, but nothing holds a height above
  # L_tip/eta_c during a run, and a trial stage can overshoot a newborn below
  # it. There the path length and the conductance go negative. The rate path
  # must refuse that as a domain error, which the stepper retries smaller,
  # whatever the leaf would have done with it: with the energy balance on, the
  # leaf refuses the negative conductance with an error the stepper cannot
  # retry, and the run ends.
  wet <- rep(0.25, 5)
  s <- TF24_Model()
  threshold <- s$pars$L_tip / ETA_C_AT_12
  for (eb in c(0, 1)) {
    s$pars$use_energy_balance <- eb
    info <- sprintf("use_energy_balance = %d", eb)
    expect_true(is.finite(
      tf24_probe(s, wet, height = threshold * (1 + 1e-6))[["assimilation"]]),
      info = info)
    for (h in c(threshold * (1 - 1e-6), threshold / 2, NaN)) {
      expect_error(tf24_probe(s, wet, height = h),
                   "^TF24: a plant of height .* shorter than its stem's terminal segment.*; it is short by",
                   class = "odelia::util::DomainError",
                   info = sprintf("%s, height = %g", info, h))
    }
  }
})

test_that("widening reaches the plant, and helps more at greater height", {
  # The end-to-end wiring check. Widening lowers resistance, so the same soil
  # supports more carbon gain -- and because the mechanism is a height exponent
  # and not a constant, the advantage must GROW with height. That ratio is the
  # thing under test; a uniform conductance rescaling would raise both equally.
  #
  # Compared at a FIXED K_s: `flat` is the shipped strategy with the profile
  # switched off, so the only difference is the height dependence. (The shipped
  # default pairs widening with a reparameterised K_s, which is a different
  # comparison -- that one is the rotation asserted above.)
  wet <- rep(0.25, 5)

  wide <- TF24_Model()
  flat <- TF24_Model()
  flat$pars$D_c <- 0
  flat$pars$L_tip <- 0

  gain <- function(h) {
    a0 <- tf24_probe(flat, wet, height = h)[["assimilation"]]
    a1 <- tf24_probe(wide, wet, height = h)[["assimilation"]]
    expect_true(is.finite(a0) && is.finite(a1))
    expect_gt(a1, a0)
    a1 / a0
  }

  expect_gt(gain(10), gain(1))
})
