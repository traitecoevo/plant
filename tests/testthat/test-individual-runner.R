context("IndividualRunner")

strategy_types <- get_list_of_strategy_types()
environment_types <- get_list_of_environment_types()

test_that("IndividualRunner", {
  for (x in names(strategy_types)) {

    e <- environment_types[[x]]
    p <- Individual(x, e)(strategy_types[[x]]())
    env <- test_environment(x, 10)
    p$compute_rates(env)

    pr <- IndividualRunner(x, e)(p, env)
    expect_is(pr, sprintf("IndividualRunner<%s,%s>",x,e))
    expect_is(pr$individual, sprintf("Individual<%s,%s>",x,e))

    expect_equal(pr$individual$internals, p$internals)

    ## This going to work with a *copy* of pr; so that won't propagate
    ## back.
    runner <- OdeRunner(x)(pr)
    expect_is(runner, "OdeRunner")
    expect_is(runner, sprintf("OdeRunner<%s>", x))
    expect_equal(runner$time, 0.0)

    expect_equal((get_individual_internals_fun(p))(runner), p$internals)
    
    continue_if <- function(obj) {
      obj$state[[1]] < 15
    }
    observer <- function(obj) {
      c(obj$time, obj$state)
    }
    pr <- IndividualRunner(x, e)(Individual(x, e)(strategy_types[[x]]()), env)
    runner <- OdeRunner(x)(pr)
    ret <- list(observer(runner))
    while (continue_if(runner)) {
      message(runner$time)
      runner$step()
      ret <- c(ret, list(observer(runner)))
    }

    m <- do.call("rbind", ret)
    colnames(m) <- c("time", runner$object$individual$ode_names)

    if (interactive()) {
      plot(height ~ time, as.data.frame(m))
    }
  }
})

test_that("get_plant_internals_fun", {
  for (x in names(strategy_types)) {
    e <- environment_types[[x]]
    p <- Individual(x, e)(strategy_types[[x]]())
    env <- test_environment(x, 10)
    p$compute_rates(env)

    runner <- OdeRunner(x)(IndividualRunner(x, e)(p, env))
    h0 <- runner$object$individual$state("height")
    runner$step()
    runner$step()
    runner$step()
    runner$step()
    h1 <- runner$object$individual$state("height")
    expect_gt(h1, h0) ## test that plants grow
  }
})

test_that("grow_individual_to_size", {
  for (x in names(strategy_types)) {
    #cat(x)
    env <- test_environment(x, 10)

    heights <- seq(1, 10)
    
    e <- environment_types[[x]]
    s <- strategy_types[[x]]()

    pp <- Individual(x, e)(s)

    if(grepl("K93", x))
      heights <- subset(heights, heights >  pp$strategy$height_0)

    res <- grow_individual_bracket(pp, heights, "height", env)

    expect_identical(res$t0, res$time[res$index])
    expect_identical(res$t1, res$time[res$index + 1L])
    expect_identical(res$y0, res$state[res$index,])
    expect_identical(res$y1, res$state[res$index + 1L,])
    ## We really do bracket the size:
    expect_true(all(res$y0[,"height"] < heights))
    expect_true(all(res$y1[,"height"] > heights))
    expect_is(res$runner, sprintf("OdeRunner<%s>", x))

    ## Then, do the search for a single case:
    i <- 3L
    tmp <- grow_individual_bisect(res$runner,
                             heights[[i]], "height",
                             res$t0[[i]], res$t1[[i]], res$y0[i,])

    ## The plant lies within the time range expected:
    expect_gte(tmp$time, res$t0[[i]])
    expect_lte(tmp$time, res$t1[[i]])
 
    j1 <- match("height", names(res$y0[i,]))
    expect_true(all(tmp$state[j1] > res$y0[i,j1]))
    expect_true(all(tmp$state[j1] < res$y1[i,j1]))

    ## test using <= & >= for fecundity & mortality, as unlike other size, fecundity could be zero
    j2 <- match(c("fecundity", "mortality"), names(res$y0[i,]))
    expect_true(all(tmp$state[j2] >= res$y0[i,j2]))
    expect_true(all(tmp$state[j2] <= res$y1[i,j2]))

    if (interactive()) {
      par(mfrow=c(3,1))
      ll <- c(res$y0[i,], res$y1[i,], tmp$state)
      y_min <- min(ll)
      y_max <- max(ll)
      plot(res$y0[i,],col="blue", type="l", ylim=c(y_min, y_max), main=x)
      points(tmp$state > res$y0[i,], col = "green", pch = 2)
      points(tmp$state < res$y1[i,], col = "pink", pch = 4)
      lines(tmp$state)
      lines(res$y1[i,], col = "red", pch = 17)
      plot(tmp$state - res$y0[i,], col = "green", pch = 2)
      plot(tmp$state - res$y1[i,], col = "pink", pch = 4)
    }
    
    ## Do all plants using the proper function:
    obj <- grow_individual_to_size(Individual(x, e)(s), heights, "height", env)
    expect_is(obj$time, "numeric")
    expect_true(all(obj$time >= res$t0))
    expect_true(all(obj$time <= res$t1))

    expect_true(all(obj$state[,j1] >= res$y0[,j1]))
    expect_true(all(obj$state[,j1] <= res$y1[,j1]))
    expect_true(all(obj$state[,j2] >= res$y0[,j2]))
    expect_true(all(obj$state[,j2] <= res$y1[,j2]))

    expect_equal(length(obj$individual), length(heights))
    expect_true(all(sapply(obj$individual, inherits, sprintf("Individual<%s,%s>",x,e))))
    expect_equal(sapply(obj$individual, function(p) p$state("height")), heights, tolerance=1e-6)
  }
})

## TODO: another useful function could be to construct splines for
## arbitrary variables during a run; we end up with all the state here
## so that should be fairly straightforward.
test_that("grow_individual_to_size", {
  for (x in names(strategy_types)) {
    strategy <- strategy_types[[x]]()
    e <- environment_types[[x]]
    pl <- Individual(x, e)(strategy)
    sizes <- c(1, 5, 10, 12, strategy$hmat)
    if(grepl("K93", x)) 
      sizes <- c(2.5, 5, 10, 12)
    env <- fixed_environment(x, 1.0)
    res <- grow_individual_to_size(pl, sizes, "height", env, 10000)

    expect_equal(res$state[, "height"], sizes, tolerance=1e-4)

    sizes2 <- c(sizes, last(sizes) * 2)
    if(x == "FF16") {
      expect_warning(res2 <- grow_individual_to_size(pl, sizes2, "height", env, 100),
                "Time exceeded time_max")
      expect_equal(length(res2$time), length(sizes2))
      expect_equal(last(res2$time), NA_real_)
      expect_false(any(is.na(res2$time[-length(sizes2)])))

      expect_silent(res3 <- grow_individual_to_size(pl, sizes2, "height", env,
                                           100, warn=FALSE))
      expect_equal(res3, res2)

      expect_silent(res4 <- grow_individual_to_size(pl, sizes2, "height", env,
                                           100, warn=FALSE, filter=TRUE))

      ## Manually filter:
      cmp <- res2
      i <- !is.na(cmp$time)
      expect_equal(res4$time, cmp$time[i])
      expect_equal(res4$plant, cmp$plant[i])
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

test_that("grow_individual_to_time", {
  for (x in names(strategy_types)) {
    strategy <- strategy_types[[x]]()
    e <- environment_types[[x]]
    pl <- Individual(x, e)(strategy)
    env <- fixed_environment(x, 1.0)
    times <- c(0, 10^(-4:3))
    res <- grow_individual_to_time(pl, times, env)
    expect_is(res$individual, "list")
    expect_equal(length(res$individual), length(times))

    expect_is(res$state, "matrix")
    expect_equal(colnames(res$state), pl$ode_names)
    expect_equal(nrow(res$state), length(times))

    expect_true(all(diff(res$state[, "height"]) > 0))

    expect_identical(res$time, times)
  }
})

test_that("Sensible behaviour on integration failure", {
  pl <- FF16_Individual()
  hyperpar <- make_FF16_hyperpar()

  env <- fixed_environment("FF16", 1)
  sizes <- seq_range(c(pl$state("height"), 50), 50)
  expect_warning(res <- grow_individual_to_size(pl, sizes, "height", env, 10, warn = TRUE, filter = TRUE),
                  "Time exceeded time_max")
  expect_is(res$individual, "list")
  expect_equal(nrow(res$state), length(res$time))
  expect_equal(res$state[,"height"], sizes[seq_len(length(res$time))], tol = 1E-5)

  ## As reported in issue #174
  traits <- trait_matrix(
    c(10.8552329005728,0.0105249489979936,633.641169104633,0.00227017846696535,3.61377429973252,0.66734,1.99575855688525,0.0567301588978768,900.481925352216,0.278101287831634,0.0183902655344452,0.523884989942541,4.34167189261234,4.520285,0.824983001217059,0.0344755338576102,0.245619588820917,0.260495,5031.78411788144,1.06992910086522,1.1368443297711,50,1),
    c("eta","lma","rho","theta","a_l1","a_l2","a_r1","a_b1","r_r","k_b","k_r","omega","B_kl1","B_kl2","B_ks1","narea","B_lf1","B_lf5","B_lf4","B_rs1","B_rb1","hmat","c_r1")
  )

  s <- strategy(traits, scm_base_parameters("FF16"), hyperpar, 1.0)
  pl <- FF16_Individual(s)

  env <- fixed_environment("FF16", 1)
  sizes <- seq_range(c(pl$state("height"), 50), 50)
  expect_warning(res <- grow_individual_to_size(pl, sizes, "height", env, 1000, warn = TRUE, filter = TRUE),
                  "50 larger sizes dropped")
  expect_equal(res$individual, list())
  expect_equal(res$time, numeric(0))
  expect_equal(nrow(res$state), 0)
})
