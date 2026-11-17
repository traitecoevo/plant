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

test_that("gd_point() gives partials, and chaining them is stand_gradient()'s total", {
  ## The guide sets the partial beside the total. A gd_point() that took the
  ## gradient through the trade-offs `p` carries would make both columns the total.
  skip_on_cran()
  source(guide_helpers("gradient_demo_helpers.R"), local = TRUE)
  tr <- unlist(gd_traits())
  pt <- gd_point(tr, lifetime = 1, refine = FALSE)
  expect_true(all(c("1.k_l", "1.r_l") %in% colnames(pt$gradient)))
  total <- gd_trait_gradient(pt$gradient, gd_hyperpar_jacobian(tr))
  metric <- "mass_above_ground"
  ## lma derives k_l and r_l, so its total differs from its partial; hmat
  ## derives nothing, so its two agree.
  expect_gt(abs(total[metric, "lma"] - pt$gradient[metric, "1.lma"]),
            1e-3 * abs(pt$gradient[metric, "1.lma"]))
  expect_equal(total[metric, "hmat"], pt$gradient[metric, "1.hmat"])
  ## The chained total is what stand_gradient() reports through the trade-offs.
  solver <- run(gd_parameters(tr, lifetime = 1), Environment("TF24"),
                Control(node_density_coordinate = "birth_date"),
                refine_schedule = FALSE, record = "states")$solver
  through <- stand_gradient(solver, metrics = metric)$gradient
  expect_equal(total[metric, "lma"], through[metric, "1.lma"], tolerance = 1e-5)
})
