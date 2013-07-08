source("helper-tree.R")

context("Load reference data")

path <- "~/Documents/Projects/veg/ebt/output_ref/T1"
## Loading the reference data takes about 1s.
ref <- tree:::load.reference.output(path)

## TODO:
##   * export core parameters
##   * export mean disturbance interval
##   * export c_d0 with more accuracy
##   * export b and c_s0

test_that("Loaded reference object seems sane", {
  expect_that(names(ref),
              equals(c("stand", "fitness", "patch_age", "params", "popn")))
  expect_that(length(ref$popn), equals(1))
  expect_that(names(ref$popn[[1]]),
              equals(c("age", "coh_m", "coh_n", "coh_r",
                       "bound_m", "bound_s", "bound_n", "bound_r",
                       "d_bound_m", "d_bound_s", "d_bound_r",
                       "popn")))
})

s <- new(Strategy)
pars <- s$parameters
ref$params <- as.list(ref$params)

test_that("Mutually missing values are as expected", {
  ## In tree, we don't include entries for:
  ##   eta_c       -- computed
  ##   a4, B5      -- not used at present
  ##   c_d4, c_d5  -- not used at present
  ##   c_ext, Pi_0 -- in parameters
  expect_that(setdiff(names(ref$params), names(pars)),
              equals(c("eta_c", "a5", "B5",
                       "c_ext", "Pi_0",
                       "c_d4", "c_d5")))

  ## We don't have:
  ##   b    -- bark area per sapwood area
  ##   c_s0 -- seedling mortality
  ##   hmat, lma, rho, s -- core traits, set elsewhere
  ##
  ## Though not exported, 'b' is set to 0.17, and 'c_s0' to 0.1 in the
  ## code (Strategy.cpp), which is the same value that we use in tree.
  expect_that(setdiff(names(pars), names(ref$params)),
              equals(c("b", "c_s0", "hmat", "lma", "rho", "s")))
})

test_that("Parameter values agree", {
  common <- intersect(names(pars), names(ref$params))

  ## The error is higher in c_d0 than in other parameters because it
  ## is exported with lower precision than it is generated with. Check
  ## the value is what we expect, then assign the new value so that
  ## the equal test works.
  expect_that(pars$c_d0, equals(ref$params$c_d0, tolerance=1e-6))
  expect_that(abs(pars$c_d0 - ref$params$c_d0),
              is_greater_than(1e-7))
  ref$params$c_d0 <- pars$c_d0
  expect_that(pars[common],
              equals(ref$params[common]))
})
