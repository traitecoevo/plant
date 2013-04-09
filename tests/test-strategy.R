source("helper-tree.R")

context("Strategy")

s <- new(Strategy)

obj <- s$get_parameters()

expected <- list(
  B1     = 0.306,
  B2     = 1.75,
  B4     = 1.71,
  Y      = 0.7,
  a1     = 5.44,
  a2     = 6.67e-05,
  a3     = 0.07, 
  a4     = 0.0286,
  b      = 0.17,
  c_Rb   = 8024,
  c_Rl   = 21000,
  c_Rr   = 217,
  c_Rs   = 4012,
  c_acc  = 4,
  c_bio  = 0.0245,
  c_d0   = 0.520393415085166, 
  c_d1   = 0.0065,
  c_d2   = 5.5,
  c_d3   = 20,
  c_p1   = 150.36,
  c_p2   = 0.19,
  c_r1   = 1,
  c_r2   = 50,
  c_s0   = 0.1,
  eta    = 12,
  hmat   = 16.5958691,
  k_b    = 0.2,
  k_r    = 1,
  lma    = 0.1978791,
  n_area = 0.00187,
  rho    = 608,
  s      = 3.8e-5,
  theta  = 4669)

keys <- sort(names(expected))
expect_that(sort(names(obj)),
            is_identical_to(sort(names(expected))))
expect_that(obj[keys], equals(expected[keys]))

## Add some new parameters:
new1 <- list(c_acc=4.1, lma=1.2)
s$set_parameters(new1)
expect_that(s$get_parameters()[names(new1)], is_identical_to(new1))
expect_that(s$get_parameters(), equals(modifyList(expected, new1)))

## Generate a failure:
new2 <- list(unknown_key=1)
expect_that(s$set_parameters(new2), throws_error())

## And have the list remain unchanged and valid
expect_that(s$get_parameters(), equals(modifyList(expected, new1)))

## TODO: error on not the first argument would not satisfy that
## property.

## Empty list should be accepted and leave things unchanged.
obj <- s$get_parameters()
s$set_parameters(list())
expect_that(s$get_parameters(), is_identical_to(obj))

## As should NULL, which is converted by R to list()
obj <- s$get_parameters()
s$set_parameters(NULL)
expect_that(s$get_parameters(), is_identical_to(obj))
