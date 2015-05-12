context("Adaptive interpolator")

test_that("Simple test case", {
  target <- function(x) sin(2*x)
  r <- c(0, 2*pi)
  s <- test_adaptive_interpolator(target, r[1], r[2])

  expect_that(s, is_a("Interpolator"))

  expect_that(s$size, equals(241))
  expect_that(nrow(s$xy), equals(s$size))

  xx_eval <- s$xy[,1]
  expect_that(s$xy[,2], is_identical_to(target(xx_eval)))

  xx_mid <- (xx_eval[-1] + xx_eval[-length(xx_eval)]) / 2
  yy_mid <- target(xx_mid)
  zz_mid <- s$eval(xx_mid)
  expect_that(zz_mid, equals(yy_mid, tolerance=2e-8))
  err <- pmax(abs(zz_mid - yy_mid), abs(1 - zz_mid / yy_mid))
  expect_that(all(err < 1e-6), is_true())
})
