
strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("run() to times grows a copy of the plant", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Individual(x, e)(strategy_types[[x]]())
    env <- Environment(x)
    p$compute_rates(env)
    h0 <- p$state("height")

    res <- run(p, times = c(0, 0, 1, 2), env = env)
    expect_inherits(res$individual[[1]], sprintf("Individual<%s,%s>", x, e))
    ## At its starting time the plant reported is the one passed in.
    expect_equal_internals(res$individual[[1]]$internals, p$internals)
    expect_identical(res$state[1, ], res$state[2, ])
    expect_true(all(diff(res$state[-1, "height"]) > 0))
    ## The plant passed in is not grown.
    expect_identical(p$state("height"), h0)
  }
  p <- FF16_Individual()
  expect_error(run(p, times = c(1, 0.5)), "times must be sorted")
  expect_error(run(p, times = c(1, NA)), "Times must be positive")
  expect_error(run(p, times = numeric(0)), "At least one time")
})

test_that("run() grows an individual to sizes", {
  for (x in names(strategy_types)) {
    env <- Environment(x)
    e <- environment_types[[x]]
    s <- strategy_types[[x]]()
    pp <- Individual(x, e)(s)

    heights <- seq(1, 10)
    if (grepl("K93", x)) {
      heights <- subset(heights, heights > pp$strategy$pars$height_0)
    }

    obj <- run(pp, sizes = heights, size_name = "height", env = env)
    expect_inherits(obj$time, "numeric")
    expect_true(all(diff(obj$time) > 0))
    expect_equal(colnames(obj$state), pp$ode_names)
    expect_equal(length(obj$individual), length(heights))
    expect_true(all(sapply(obj$individual, inherits, sprintf("Individual<%s,%s>", x, e))))

    ## Each size's time is located on the interpolant of the step that
    ## crosses it, cubic Hermite under Cash-Karp, so the height there misses
    ## the size by the rate times that interpolant's error in time: up to
    ## ~8e-6 m in FF16. The largest miss is checked, not a pooled mean.
    expect_lt(max(abs(obj$state[, "height"] - heights)), 1e-5)

    ## Time, state and individual describe the same instant: the individual
    ## holds the state, and a run to those times, landing on each in turn,
    ## holds it too, to the bit.
    for (i in seq_along(heights)) {
      expect_identical(obj$individual[[i]]$ode_state, unname(obj$state[i, ]))
    }
    expect_identical(run(pp, times = obj$time, env = env)$state, obj$state)

    ## The trajectory starts where the plant did and holds each size's step.
    expect_equal(colnames(obj$trajectory), c("time", pp$ode_names))
    expect_identical(unname(obj$trajectory[1, ]), c(0, pp$ode_state))
    expect_true(all(diff(obj$trajectory[, "time"]) > 0))
    k <- match(obj$time, obj$trajectory[, "time"])
    expect_false(anyNA(k))
    expect_identical(unname(obj$trajectory[k, -1]), unname(obj$state))
  }
})

test_that("run() to sizes reports a repeated size, or one already reached, at once", {
  x <- "FF16"
  env <- Environment(x)
  pp <- Individual(x, "FF16_Env")(FF16_Model())
  h0 <- pp$state("height")
  obj <- run(pp, sizes = c(h0, 2, 2, 3), env = env)
  expect_identical(obj$time[[1]], 0)
  expect_identical(obj$time[[2]], obj$time[[3]])
  expect_identical(obj$state[2, ], obj$state[3, ])
  expect_gt(obj$time[[4]], obj$time[[3]])

  ## A landing can stop a hair past its size; a next size inside that hair
  ## is reached where the run stands.
  d <- run(pp, sizes = 5, env = env)$state[1, "height"] - 5
  if (d > 0) {
    obj <- run(pp, sizes = c(5, 5 + d / 2), env = env)
    expect_identical(obj$time[[2]], obj$time[[1]])
  }

  expect_error(run(pp, sizes = 2, size_name = "leaf_area", env = env),
               "'leaf_area' is not one of the ODE states")
  expect_error(run(pp, sizes = c(3, 2), env = env), "sorted")
  expect_error(run(pp, sizes = c(2, Inf), env = env), "finite")
  expect_error(run(pp, sizes = h0 / 2, env = env), "already bigger")
  for (bad in list(NA_real_, NaN, -Inf, -1, c(1, 2))) {
    expect_error(run(pp, sizes = 2, env = env, time_max = bad),
                 "time_max must be a single non-negative number")
  }
})

test_that("a run to sizes that fails warns with the error and what was dropped", {
  pp <- Individual("FF16", "FF16_Env")(FF16_Model())
  pp$compute_rates(Environment("FF16"))
  ## This plant's mortality rate is not zero, so the failure is not read as
  ## possibly something other than its death, and only the first warning is
  ## given.
  expect_warning(integration_failed(pp, "the step failed", 3L, TRUE),
                 "integration failed with error: \nthe step failed\n  3 larger sizes dropped",
                 fixed = TRUE)
  expect_silent(integration_failed(pp, "the step failed", 3L, FALSE))
})

test_that("run() grows an individual to sizes, and stops at time_max", {
  for (x in names(strategy_types)) {
    strategy <- strategy_types[[x]]()
    e <- environment_types[[x]]
    pl <- Individual(x, e)(strategy)
    sizes <- c(1, 5, 10, 12, strategy$pars$hmat)
    
    if(grepl("K93", x)) 
      sizes <- c(2.5, 5, 10, 12)

    if(grepl("TF24", x)) 
      sizes <- c(1, 5, 10, 12)
    env <- Environment(x)
    env$set_fixed_environment(1.0, height_max = 150)
    
    res <- run(pl, sizes = sizes, size_name = "height", env = env, time_max = 10000)

    expect_equal(res$state[, "height"], sizes, tolerance=1e-4)

    sizes2 <- c(sizes, dplyr::last(sizes) * 2)
    if(x %in% c("FF16", "TF24")){
      expect_warning(res2 <- run(pl, sizes = sizes2, size_name = "height", env = env,
                                 time_max = 100),
                "Time reached time_max")
      expect_equal(length(res2$time), length(sizes2))
      expect_equal(dplyr::last(res2$time), NA_real_)
      expect_false(any(is.na(res2$time[-length(sizes2)])))

      expect_silent(res3 <- run(pl, sizes = sizes2, size_name = "height", env = env,
                                time_max = 100, warn=FALSE))
      ## warn=FALSE must produce the same numeric result as the warning version.
      ## Compare the value payload only: the `individual`/`env` entries are live
      ## reference objects whose C++ `.ptr` always differs between two calls
      ## (edition-3 waldo would flag those pointers; edition 2 ignored them).
      expect_equal(res3$time, res2$time)
      expect_equal(res3$state, res2$state)
      expect_equal(res3$trajectory, res2$trajectory)

      expect_silent(res4 <- run(pl, sizes = sizes2, size_name = "height", env = env,
                                time_max = 100, warn=FALSE, filter=TRUE))

      ## Manually filter:
      cmp <- res2
      i <- !is.na(cmp$time)
      expect_equal(res4$time, cmp$time[i])
      expect_equal(lapply(res4$individual, function(p) p$ode_state),
                   lapply(cmp$individual[i], function(p) p$ode_state))
      expect_equal(res4$state, cmp$state[i,])
      expect_equal(res4$trajectory, cmp$trajectory)

      if (FALSE) {
        plot(height ~ time, as.data.frame(res$trajectory), type="l")
        points(res$time, res$state[, "height"], pch=19)

        plot(height ~ time, as.data.frame(res2$trajectory), type="l")
        points(res2$time, res2$state[, "height"], pch=19)
      }
    }
  }
})

test_that("run() grows an individual to times", {
  for (x in names(strategy_types)) {
    strategy <- strategy_types[[x]]()
    e <- environment_types[[x]]
    pl <- Individual(x, e)(strategy)
    env <- Environment(x)
    env$set_fixed_environment(1.0, height_max = 150)

    times <- c(0, 10^(-4:3))
    res <- run(pl, times = times, env = env)
    expect_inherits(res$individual, "list")
    expect_equal(length(res$individual), length(times))

    expect_inherits(res$state, "matrix")
    expect_equal(colnames(res$state), pl$ode_names)
    expect_equal(nrow(res$state), length(times))

    expect_true(all(diff(res$state[, "height"]) > 0))

    expect_identical(res$time, times)
  }
})

test_that("Sensible behaviour on integration failure", {
  pl <- FF16_Individual()
  hyperpar <- default_hyperpar("FF16")

  env <- Environment("FF16")
  env$set_fixed_environment(1.0, height_max = 150)

  sizes <- seq_range(c(pl$state("height"), 50), 50)
  expect_warning(res <- run(pl, sizes = sizes, size_name = "height", env = env,
                            time_max = 10, warn = TRUE, filter = TRUE),
                  "Time reached time_max")
  expect_inherits(res$individual, "list")
  expect_equal(nrow(res$state), length(res$time))
  expect_lt(max(abs(res$state[, "height"] - sizes[seq_along(res$time)])), 1e-5)

  ## As reported in issue #174
  ##
  ## ⚠️ THIS MATRIX USED TO CARRY TEN MORE COLUMNS, NINE OF WHICH DID NOTHING:
  ## `B_kl1`, `B_kl2`, `B_ks1`, `B_lf1`, `B_lf4`, `B_lf5`, `B_rs1`, `B_rb1` are
  ## CONSTRUCTOR arguments of the FF16 trade-off function, not traits -- the hyperpar is
  ## built with no arguments on the line above, so naming them here never
  ## configured anything -- and `c_r1` is not a parameter of this model at all.
  ## They were appended to `pars` as junk elements and ignored by the generated
  ## Rcpp converter. #636 refuses an unknown trait name, which is what surfaced
  ## them.
  ##
  ## `narea` is KEPT: FF16's hyperpar genuinely reads it from the matrix, to
  ## derive `a_p1`, `a_p2` and `r_l`. Measured before removing the rest -- of the
  ## ten, `narea` alone moved any parameter of the strategy -- so this run is
  ## unchanged and the #174 regression below still exercises what it did.
  traits <- trait_matrix(
    c(10.8552329005728,0.0105249489979936,633.641169104633,0.00227017846696535,3.61377429973252,0.66734,1.99575855688525,0.0567301588978768,900.481925352216,0.278101287831634,0.0183902655344452,0.523884989942541,0.0344755338576102,50),
    c("eta","lma","rho","theta","a_l1","a_l2","a_r1","a_b1","r_r","k_b","k_r","omega","narea","hmat")
  )

  s <- generate_strategy(testing_parameters("FF16"), traits, hyperpar = hyperpar, birth_rate = 1.0)[[1]]
  pl <- FF16_Individual(s)

  env <-Environment("FF16")
  env$set_fixed_environment(1.0, height_max = 150)
     
  sizes <- seq_range(c(pl$state("height"), 50), 50)
  ## This plant never grows. The first size is its starting height, reached
  ## at the start; it reaches none of the others.
  expect_warning(res <- run(pl, sizes = sizes, size_name = "height", env = env,
                            time_max = 1000, warn = TRUE, filter = TRUE),
                  "Time reached time_max, 49 larger sizes dropped")
  expect_equal(length(res$individual), 1L)
  expect_identical(res$time, 0)
  expect_equal(nrow(res$state), 1)
})
