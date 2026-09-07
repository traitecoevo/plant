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
# plant/stem_hydraulics.h and notes/plan-tf24-height-hydraulics.md.
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
  p <- TF24_Strategy()$pars
  # Widening is on; theta_c stays at zero so `theta` keeps its whole-plant
  # meaning and no parameter-file migration is needed yet.
  expect_identical(p$D_c, 0.2)
  expect_identical(p$theta_c, 0)
  expect_identical(p$L_tip, 0.02)
  # K_s is now a TERMINAL-SEGMENT conductivity, back-derived from the old
  # whole-stem 1 so that resistance is unchanged at the anchor height. Asserted
  # against the helper rather than a literal, so the two cannot drift apart.
  # Deliberately expect_equal, not expect_identical. The header default is a
  # literal while the helper computes it through log() and expm1(), whose
  # last-ULP behaviour is libm-dependent -- so a bitwise comparison here would be
  # asserting something about the platform's maths library rather than about
  # plant. 1e-12 is far tighter than any real drift and still catches the thing
  # that matters, which is the literal and the helper falling out of step.
  expect_equal(p$K_s, TF24_K_s_from_whole_stem(1), tolerance = 1e-12)
  expect_equal(1 / p$K_s, 2.9781957045054700, tolerance = 1e-12)
})

test_that("zeroing the profile parameters recovers the height-linear model", {
  # Invariance criterion I7 survives the reparameterisation as an explicit
  # rather than as the default. expect_identical, not expect_equal: a value that
  # drifted to 1e-16 would still pass expect_equal while silently taking the
  # beta == 0 branch out of play, and the SCM pins elsewhere sit at tolerance
  # 2e-2, which would absorb the consequence without anyone noticing.
  s <- TF24_Strategy()
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
    expect_identical(s$pars$K_s * s$pars$theta / L,
                     s$pars$K_s * s$pars$theta / (h * eta_c))
  }
})

test_that("resistance is unchanged at the anchor height and rotates about it", {
  # The single point of agreement between the old and new models. R_L depends on
  # theta and K_s only through their ratio, so there is exactly one free scalar
  # and agreement holds at exactly ONE height -- never two. Below the anchor
  # every plant is more resistant than under the old model, above it every plant
  # is less: the reparameterisation is a rotation, not a rescaling.
  p <- TF24_Strategy()$pars
  eta_c <- 1 - 2 / (1 + p$eta) + 1 / (1 + 2 * p$eta)
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
  s <- TF24_Strategy()
  m <- trait_matrix(s$pars$lma, "lma")
  derived <- TF24_hyperpar(m, s, filter = FALSE)
  expect_equal(derived[, "stem_P50"][[1]], 2.888725665336019, tolerance = 1e-10)

  # And the un-re-anchored value, recorded so the size of the averted leak stays
  # visible: this is what stem_P50 would be at B_Hv1 = 0.4607063.
  expect_equal(10^(0.4607063 - 0.2 * log10(s$pars$K_s)), 3.5933287036419,
               tolerance = 1e-10)
})

test_that("stem path parameters are settable", {
  s <- TF24_Strategy()
  s$pars$D_c <- 0.3
  s$pars$L_tip <- 0.05
  expect_equal(s$pars$D_c, 0.3)
  expect_equal(s$pars$L_tip, 0.05)
})

test_that("a non-zero theta_c is refused rather than half-applied", {
  # theta is not a hydraulics-only trait: it also sets area_sapwood, area_bark,
  # mass_sapwood (hence construction cost, respiration, turnover and NSC
  # capacity) and the hard-coded dmass_sapwood_darea_leaf derivative, all of
  # which read a flat pars.theta. Profiling it on the hydraulic side alone would
  # give a plant that conducts as though theta varied along the stem and is
  # built as though it did not -- two different plants sharing one trait.
  #
  # The field is declared so the closed form can be exercised and swept in
  # isolation, but a strategy carrying it must not build. Invariance criterion
  # I13. Delete this test in the same change that profiles theta everywhere.
  s <- TF24_Strategy()
  s$pars$theta_c <- 0.4
  expect_error(TF24_Individual(s), "theta_c is not implemented")

  # ...including when it would exactly cancel the widening term, which is the
  # case a beta-only guard would miss.
  s2 <- TF24_Strategy()
  s2$pars$D_c <- 0.2
  s2$pars$theta_c <- -0.4
  expect_identical(2 * s2$pars$D_c + s2$pars$theta_c, 0)
  expect_error(TF24_Individual(s2), "theta_c is not implemented")

  # Zero is fine, and widening on its own is unaffected -- it enters through
  # k_s(L), which has no structural counterpart to keep in step.
  expect_silent(TF24_Individual(TF24_Strategy()))
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

test_that("theta_c enters with the sign that compensates", {
  # theta(L) = theta * (L/L_tip)^(-theta_c): theta FALLS basipetally while the
  # Huber value 1/theta rises. A positive theta_c must therefore reduce the
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
  s <- TF24_Strategy()
  s$pars$L_tip <- 0
  expect_error(TF24_Individual(s), "L_tip")

  s$pars$L_tip <- 0.02
  expect_silent(TF24_Individual(s))

  # Longer than the birth-size flow path: the plant would be shorter than one
  # terminal segment, the path length would go negative, and the symptom
  # without this guard is an unattributable NaN inside the leaf solver.
  s$pars$L_tip <- 100
  expect_error(TF24_Individual(s), "terminal segment")

  # The guard is scoped to an active profile: L_tip is ignored at beta == 0.
  inert <- TF24_Strategy()
  inert$pars$D_c <- 0
  inert$pars$L_tip <- 100
  expect_silent(TF24_Individual(inert))
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

  wide <- TF24_Strategy()
  flat <- TF24_Strategy()
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
