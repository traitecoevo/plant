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

  # Possibly update to be relevant?

  # b_0 <- c(0.059)
  #ret <- K93_hyperpar(trait_matrix(b_0, "b_0"), s)

 #s expect_true(all(c("b_0", "b_1", "b_2") %in% colnames(ret)))
  #expect_equal(ret[, "b_0"], b_0)
  #expect_equal(ret[, "b_0"], c(1.46678,0.028600), tolerance=1e-5)

  ## Empty trait matrix:
  ret <- K93_hyperpar(trait_matrix(numeric(0), "b_0"), s)
  expect_equal(ret, trait_matrix(numeric(0), "b_0"))
})

