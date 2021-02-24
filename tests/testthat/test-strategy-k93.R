# Built from  tests/testthat/test-strategy-ff16r.R on Wed Aug 12 15:33:08 2020 using the scaffolder, from the strategy:  FF16r
# Built from  tests/testthat/test-strategy-ff16.R on Fri Jul  3 08:14:35 2020 using the scaffolder, from the strategy:  FF16
context("Strategy-K93")

test_that("Defaults", {
  expected <- list(
   height_0 = 2.0,
   b_0 = 0.059,
   b_1 = 0.012,
   b_2 = 0.00041,
   c_0 = 0.008,
   c_1 = 0.00044,
   d_0 = 0.00073,
   d_1 = 0.044,
   S_D = 1,
   control = Control())

  keys <- sort(names(expected))

  s <- K93_Strategy()
  expect_is(s, "K93_Strategy")

  expect_identical(sort(names(s)), keys)
  expect_identical(unclass(s)[keys], expected[keys])
})

test_that("K93 collect_all_auxillary option", {

  s <- K93_Strategy()
  p <- K93_Individual(s)
  expect_equal(p$aux_size, 1)
  expect_equal(length(p$internals$auxs), 1)
  expect_equal(p$aux_names, c(
    "competition_effect"
  ))
})

test_that("Reference comparison", {
  s <- K93_Strategy()
  p <- K93_Individual(s)

  expect_identical(p$strategy, s)

  ## Set the height to something (here 10)
  h0 <- 10
  p$set_state("height", h0)

  expect_identical(p$state("height"), h0)

  ## Check: Is this redundant now
  ## We now use
  vars <- p$internals
  expect_identical(p$state("height"), vars$states[which(p$ode_names == "height")])
})

test_that("Critical Names", {
  s <- K93_Strategy()
  my_names <- K93_Individual(s)$ode_names
  expect_identical(my_names[1:3], c("height", "mortality", "fecundity"))
})

test_that("K93_Strategy hyper-parameterisation", {
  s <- K93_Strategy()

  ## Hyperpars should just pass through:
  ret <- K93_hyperpar(trait_matrix(numeric(0), "b_0"), s)
  expect_equal(ret, trait_matrix(numeric(0), "b_0"))
})

## Number of ODE steps is unstable - needs review
# test_that("K93 seed rain is unchanged", {

#  # Generic parameters
#  p0 <- scm_base_parameters("K93")
#  p0$k_I <- 1e-6
#  p0$disturbance_mean_interval <- 10

#  # Use single sp. defaults
#  p1 <- expand_parameters(trait_matrix(0.059, "b_0"), p0, mutant = FALSE)
#  p1$seed_rain <- 20

#  out <- run_scm(p1)
#  expect_equal(out$seed_rains, 0.0752, tolerance = 1e-4)
#  expect_equal(out$ode_times[c(10, 100)], c(0.000070, 4.500004), tolerance = 1e-5)

#  # Three species from paper
#  sp <- trait_matrix(c(0.042, 0.063, 0.052,
#                       8.5e-3, 0.014, 0.015,
#                       2.2e-4, 4.6e-4, 3e-4,
#                       0.008, 0.008, 0.008,
#                       1.8e-4, 4.4e-4, 5.1e-4,
#                       1.4e-4, 2.5e-3, 8.8e-3,
#                       0.044, 0.044, 0.044),
#                     c("b_0", "b_1", "b_2",
#                       "c_0", "c_1", "d_0", "d_1"))

#  p2 <- expand_parameters(sp, p0, mutant = FALSE)
#  p2$seed_rain <- c(20, 20, 20)
#  out <- run_scm(p2)

#  expect_equal(out$seed_rains, c(0.0025, 0.2314, 0.2195), tolerance = 1e-4)
#  expect_equal(length(out$ode_times), 224)
#})
