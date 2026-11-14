test_that("gd_point() reads a value for every gradient column", {
  ## The gradient has a column for the constant birth rate, which is not a
  ## parameter: its value is the strategy's birth_rate_y.
  skip_on_cran()
  source(guide_helpers("gradient_demo_helpers.R"), local = TRUE)
  pt <- gd_point(unlist(gd_traits()), lifetime = 1, refine = FALSE)
  expect_identical(names(pt$x), colnames(pt$gradient))
  expect_true(all(is.finite(pt$x)))
  p <- gd_parameters(unlist(gd_traits()), lifetime = 1)
  expect_identical(pt$x[["1.birth_rate"]], p$strategies[[1]]$birth_rate_y)
  ## Every other column is a parameter, read from the strategy's pars.
  pars <- p$strategies[[1]]$pars
  other <- setdiff(names(pt$x), "1.birth_rate")
  expect_identical(unname(pt$x[other]),
                   vapply(sub("^1\\.", "", other), function(nm) as.numeric(pars[[nm]]),
                          numeric(1), USE.NAMES = FALSE))
})
