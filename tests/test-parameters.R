source("helper-tree.R")

context("Parameters")

p <- new(Parameters)

obj <- p$get_params()

expected <- list(
  B1     = 0.306,
  B2     = 1.75,
  B4     = 1.71,
  Pi_0   = 0.25,
  Y      = 0.7,
  a1     = 5.44,
  a2     = 6.67e-05,
  a3     = 0.07, 
  a4     = 0.0286,
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
  eta    = 12,
  hmat   = NA_real_,
  k_b    = 0.2,
  k_r    = 1,
  lma    = NA_real_,
  n_area = 0.00187,
  rho    = NA_real_,
  s      = NA_real_,
  theta  = 4669)
core <- c("lma", "rho", "hmat", "s")

keys <- sort(names(expected))
expect_that(sort(names(obj)),
            is_identical_to(sort(names(expected))))
expect_that(obj[keys], equals(expected[keys]))
expect_that(all(sapply(obj[core], is.na)), is_true())
expect_that(all(!sapply(obj[setdiff(keys, core)], is.na)), is_true())
