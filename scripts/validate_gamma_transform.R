# Issue #468 validation: closed-form Weibull hydraulic transform vs the
# numerical implementation used in the code.
#
# The supply-side transpiration integrates the xylem vulnerability curve
#   K(m)/Kmax = proportion_of_conductivity(m) = exp(-(m/stem_b)^stem_c),  m = -psi (MPa, >0)
# giving the hydraulic transform
#   F(m) = int_0^m exp(-(s/stem_b)^stem_c) ds.
#
# Substituting u = (s/stem_b)^stem_c gives the closed form (lower incomplete gamma):
#   F(m) = (stem_b/stem_c) * gamma_lower(1/stem_c, (m/stem_b)^stem_c)
#        = (stem_b/stem_c) * Gamma(1/stem_c) * pgamma((m/stem_b)^stem_c, 1/stem_c)        [R form]
# This matches issue #468's P50 form because the code sets
#   stem_b = stem_P50 / (ln 2)^(1/stem_c)  =>  ln2*(-psi/stem_P50)^stem_c = (m/stem_b)^c.
# (The issue's upper-incomplete-gamma form differs only by the additive constant
#  (stem_b/stem_c)Gamma(1/stem_c), which cancels in F(psi2) - F(psi1).)
#
# Compares the closed form against the two numerical paths actually in the code:
#   transpiration()                  -> production cubic-spline F (tk::spline)
#   transpiration_full_integration() -> QAG adaptive quadrature (tol 1e-3)
# and validates the inverse spline psi_from_transpiration.
#
# Run from the package root:  Rscript scripts/validate_gamma_transform.R

suppressMessages(pkgload::load_all(".", recompile = FALSE, quiet = TRUE))

f_analytic <- function(m, stem_b, stem_c) (stem_b / stem_c) * gamma(1 / stem_c) * pgamma((m / stem_b)^stem_c, 1 / stem_c)

make_leaf <- function(stem_c, stem_b, ncontrol = 100) {
  Leaf(vcmax_25 = 100, jmax_25 = 167 * 100, stem_c = stem_c, stem_b = stem_b,
       psi_crit = 5, root_c = 2.65, root_b = 1.29,
       root_psi_crit = 1.29 * (log(1 / 0.05))^(1 / 2.65),
      TF24_beta2 = 1, a = 0.3, curv_fact_elec_trans = 0.7,
       curv_fact_colim = 0.99, GSS_tol_abs = 1e-8,
       vulnerability_curve_ncontrol = ncontrol, ci_abs_tol = 1e-6,
       ci_niter = 1000, TF24_cost_scale = 46.33)
}

report <- function(label, stem_c, stem_b, ncontrol = 100) {
  l <- make_leaf(stem_c, stem_b, ncontrol)
  l$leaf_specific_conductance_max_ <- 1  # so transpiration() returns F directly

  psi_max <- stem_b * (log(1 / 0.01))^(1 / stem_c)   # spline upper bound (1% conductivity)
  # Replicate the C++ knot loop exactly: floating-point `psi += step`
  # accumulation determines the true last knot, so we stay inside the domain
  # (the spline has extrapolate = FALSE).
  step <- psi_max / ncontrol
  knots <- 0
  psi_spline <- step
  while (psi_spline <= psi_max) {
    knots <- c(knots, psi_spline)
    psi_spline <- psi_spline + step
  }
  last_knot <- knots[length(knots)]
  m <- seq(step / 10, last_knot * (1 - 1e-9), length.out = 2001)

  fa   <- f_analytic(m, stem_b, stem_c)
  fsp  <- vapply(m, function(x) l$transpiration(x, 0), numeric(1))                  # spline
  fqag <- vapply(m, function(x) l$transpiration_full_integration(x, 0), numeric(1)) # QAG

  rel <- function(a, e) max(abs(a - e) / pmax(abs(e), 1e-12))

  # Decompose spline error: AT the knots the spline reproduces the trapezoid
  # cumulative sum exactly, so knot error == trapezoid quadrature bias; off-knot
  # error adds cubic interpolation on top.
  kint  <- knots[knots > 0 & knots < last_knot]
  fk_sp <- vapply(kint, function(x) l$transpiration(x, 0), numeric(1))
  trap_bias <- max(abs(fk_sp - f_analytic(kint, stem_b, stem_c)))

  # Inverse spline: feed exact F(m), expect to recover m. Restrict to the
  # interior where E stays within the spline's knot range.
  inv_ok <- fa <= max(fsp)
  m_rec  <- vapply(fa[inv_ok], function(E) l$transpiration_to_psi_stem(E, 0), numeric(1))
  m_inv  <- m[inv_ok]

  cat(sprintf("\n=== %s  (stem_c=%.4f, stem_b=%.4f, psi_max=%.3f, ncontrol=%g) ===\n",
              label, stem_c, stem_b, psi_max, ncontrol))
  cat(sprintf("  QAG quadrature   vs analytic:  max|abs|=%.3e  max|rel|=%.3e\n",
              max(abs(fqag - fa)), rel(fqag, fa)))
  cat(sprintf("  spline F         vs analytic:  max|abs|=%.3e  max|rel|=%.3e\n",
              max(abs(fsp - fa)),  rel(fsp, fa)))
  cat(sprintf("    of which trapezoid bias (at knots)=%.3e  cubic interp adds=%.3e\n",
              trap_bias, max(abs(fsp - fa)) - trap_bias))
  cat(sprintf("  inverse spline   vs identity:  max|abs|=%.3e MPa  max|rel|=%.3e\n",
              max(abs(m_rec - m_inv)), rel(m_rec, m_inv)))
}

# Production TF24 default: stem_P50 = 1.85
p50   <- 1.85
c_def <- log(log(1 - 0.5) / log(1 - 0.88)) / (log(p50) - log(5.16))
b_def <- p50 / (-log(1 - 0.5))^(1 / c_def)
report("TF24 production default", c_def, b_def, 100)

# Leaf-test params (stem_c = 2.04, stem_b = 3): a steeper curve
report("test-leaf params", 2.04, 3.0, 100)

# Effect of spline resolution on the production default
report("production default, ncontrol=1000", c_def, b_def, 1000)
report("production default, ncontrol=30",   c_def, b_def, 30)
