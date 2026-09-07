# Extracted from test-strategy-tf24.R:408

# test -------------------------------------------------------------------------
linear <- function(p) {
    # Every strategy, not just the first: resetting only strategies[[1]] leaves
    # the second species on the path integral and silently compares two
    # different models.
    for (i in seq_along(p$strategies)) {
      s <- p$strategies[[i]]
      s$pars$D_c <- 0
      s$pars$theta_c <- 0
      s$pars$L_tip <- 0
      s$pars$K_s <- 1
      p$strategies[[i]] <- s
    }
    p
  }
p0 <- scm_base_parameters("TF24")
env <- Environment("TF24")
ctrl <- Control()
p0$max_patch_lifetime <- 5
p1 <- add_strategies(p0, trait_matrix(c(0.0825, 5), c("lma", "hmat")),
                       hyperpar = TF24_hyperpar, birth_rate = list(20))
out <- run_scm(linear(p1), env, ctrl)
expect_equal(out$offspring_production, 30.22207354, tolerance = 2e-2)
p2 <- add_strategies(p0, trait_matrix(c(0.0825, 0.10, 5, 5), c("lma", "hmat")),
                       hyperpar = TF24_hyperpar, birth_rate = list(20, 20))
out2 <- run_scm(linear(p2), env, ctrl)
expect_equal(out2$offspring_production[[1]], 23.20349831, tolerance = 2e-2)
