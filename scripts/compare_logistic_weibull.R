## Issue #467: can a logistic vulnerability curve replace the Weibull?
##
## Pure-R, no plant build required: both the curve and its forward transform F
## are closed form. Magnitude convention m = -psi >= 0, matching the +ve `psi`
## argument of Leaf::proportion_of_conductivity.
##
##   Weibull   f(m) = exp(-(m/b)^c)              [current model]
##   Logistic  g(m) = 1 / (1 + exp(a*(m - m50))) [proposed, issue #467]
##
## We ask: how closely can a logistic match plant's production Weibull curves
## (stem c~1.09, root c~2.68), how much does the *transpiration integral*
## F(m)=int_0^m f ds shift (that integral, not f itself, drives the model), and
## where do the tails diverge.

options(digits = 6)

## ---- plant production parameters (tf24_strategy.h) ----
p_50 <- 1.85
c    <- log(log(1 - 0.5) / log(1 - 0.88)) / (log(p_50) - log(5.16))
b    <- p_50 / (-log(1 - 0.5))^(1 / c)
psi_crit <- b * (log(1 / 0.05))^(1 / c)

root_c <- 2.680147
root_b <- 3.898245
root_psi_crit <- root_b * (log(1 / 0.05))^(1 / root_c)
root_p50 <- root_b * (-log(0.5))^(1 / root_c)   # magnitude at 50% loss

cat(sprintf("STEM : c=%.4f b=%.4f p50=%.4f psi_crit=%.4f\n", c, b, p_50, psi_crit))
cat(sprintf("ROOT : c=%.4f b=%.4f p50=%.4f psi_crit=%.4f\n",
            root_c, root_b, root_p50, root_psi_crit))

## ---- curves ----
weib <- function(m, b, c) exp(-(m / b)^c)
logi <- function(m, a, m50) 1 / (1 + exp(a * (m - m50)))

## Forward transform F(m) = int_0^m f(s) ds.
## Weibull: lower incomplete gamma. Logistic: elementary (the issue's selling point).
F_weib <- function(m, b, c) (b / c) * gamma(1 / c) * pgamma((m / b)^c, 1 / c)
F_logi <- function(m, a, m50)
  m - (1 / a) * log1p(exp(a * (m - m50))) + (1 / a) * log1p(exp(-a * m50))

## ---- two principled ways to choose the logistic shape `a` (m50 = p50) ----
## (1) slope-match at the 50% point: g'(m50) = -a/4, f'(p50) = -0.5*c*ln2/p50
a_slope <- function(c, p50) 2 * c * log(2) / p50
## (2) least-squares best fit of (a, m50) over the working domain [0, crit]
fit_ls <- function(b, c, crit, p50) {
  m <- seq(0, crit, length.out = 2000); y <- weib(m, b, c)
  o <- optim(c(a_slope(c, p50), p50),
             function(p) sum((logi(m, p[1], p[2]) - y)^2),
             method = "Nelder-Mead", control = list(reltol = 1e-12))
  list(a = o$par[1], m50 = o$par[2])
}

report <- function(label, b, c, crit, p50) {
  m  <- seq(0, crit, length.out = 5000); yw <- weib(m, b, c)
  Fw <- F_weib(m, b, c); Fmax <- F_weib(crit, b, c)

  a_s <- a_slope(c, p50); yl_s <- logi(m, a_s, p50)
  ls  <- fit_ls(b, c, crit, p50); yl_f <- logi(m, ls$a, ls$m50)

  dev    <- function(yl) c(max = max(abs(yl - yw)), rms = sqrt(mean((yl - yw)^2)))
  intdev <- function(a, m50) max(abs(F_logi(m, a, m50) - Fw)) / Fmax

  ds <- dev(yl_s); df <- dev(yl_f)
  cat(sprintf("\n=== %s  (domain [0, %.3f]) ===\n", label, crit))
  cat(sprintf("  slope-match : a=%.4f m50=%.4f | cond max-dev=%.4f rms=%.4f | F rel-dev=%.4f\n",
              a_s, p50, ds["max"], ds["rms"], intdev(a_s, p50)))
  cat(sprintf("  best-fit LS : a=%.4f m50=%.4f | cond max-dev=%.4f rms=%.4f | F rel-dev=%.4f\n",
              ls$a, ls$m50, df["max"], df["rms"], intdev(ls$a, ls$m50)))
}

report("STEM", b, c, psi_crit, p_50)
report("ROOT", root_b, root_c, root_psi_crit, root_p50)

cat("\n--- tail comparison: proportion of conductivity (stem, slope-matched) ---\n")
m_tail <- c(0.5, 1, 2, 3, 4, 5, 6, 7)
print(data.frame(m = m_tail,
                  weibull = weib(m_tail, b, c),
                  logistic = logi(m_tail, a_slope(c, p_50), p_50)), digits = 4)
