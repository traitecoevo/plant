context("PlantRunner")

test_that("PlantRunner", {
  p <- FFW16_PlantPlus(FFW16_Strategy())
  env <- test_environment(10)
  p$compute_vars_phys(env)

  pr <- PlantRunner(p, env)
  expect_that(pr, is_a("PlantRunner"))
  expect_that(pr$plant, is_a("FFW16_PlantPlus"))
  expect_that(pr$plant$internals, is_identical_to(p$internals))

  ## This going to work with a *copy* of pr; so that won't propagate
  ## back.
  runner <- OdeRunner("PlantRunner")(pr)
  expect_that(runner, is_a("OdeRunner"))
  expect_that(runner, is_a("OdeRunner<PlantRunner>"))
  expect_that(runner$time, equals(0.0))

  expect_that(oderunner_plant_size(runner),
              is_identical_to(p$internals))

  continue_if <- function(obj) {
    obj$state[[1]] < 15
  }
  observer <- function(obj) {
    c(obj$time, obj$state)
  }
  pr <- PlantRunner(FFW16_PlantPlus(FFW16_Strategy()), env)
  runner <- OdeRunner("PlantRunner")(pr)
  ret <- list(observer(runner))
  while (continue_if(runner)) {
    message(runner$time)
    runner$step()
    ret <- c(ret, list(observer(runner)))
  }

  m <- do.call("rbind", ret)
  colnames(m) <- c("time", runner$object$plant$ode_names)

  if (interactive()) {
    plot(height ~ time, as.data.frame(m))
  }
})

test_that("grow_plant_to_size", {
  env <- test_environment(10)
  heights <- seq(1, 10)
  s <- FFW16_Strategy()

  res <- grow_plant_bracket(FFW16_PlantPlus(s), heights, "height", env)

  expect_that(res$t0, is_identical_to(res$time[res$index]))
  expect_that(res$t1, is_identical_to(res$time[res$index + 1L]))
  expect_that(res$y0, is_identical_to(res$state[res$index,]))
  expect_that(res$y1, is_identical_to(res$state[res$index + 1L,]))
  ## We really do bracket the size:
  expect_that(all(res$y0[,"height"] < heights), is_true())
  expect_that(all(res$y1[,"height"] > heights), is_true())
  expect_that(res$runner, is_a("OdeRunner<PlantRunner>"))

  ## Then, do the search for a single case:
  i <- 3L
  tmp <- grow_plant_bisect(res$runner,
                           heights[[i]], "height",
                           res$t0[[i]], res$t1[[i]], res$y0[i,])

  ## The plant lies within the range expected:
  expect_that(tmp$time, is_within_interval(res$t0[[i]], res$t1[[i]]))
  expect_that(all(tmp$state > res$y0[i,]), is_true())
  expect_that(all(tmp$state < res$y1[i,]), is_true())

  ## Do all plants using the proper function:
  obj <- grow_plant_to_size(FFW16_PlantPlus(s), heights, "height", env)
  expect_that(obj$time, is_a("numeric"))
  expect_that(all(obj$time > res$t0), is_true())
  expect_that(all(obj$time < res$t1), is_true())

  expect_that(all(obj$state > res$y0), is_true())
  expect_that(all(obj$state < res$y1), is_true())

  expect_that(length(obj$plant), equals(length(heights)))
  expect_that(all(sapply(obj$plant, inherits, "FFW16_PlantPlus")), is_true())
  expect_that(sapply(obj$plant, function(p) p$height),
              equals(heights, tolerance=1e-6))
})

## TODO: another useful function could be to construct splines for
## arbitrary variables during a run; we end up with all the state here
## so that should be fairly straightforward.

test_that("grow_plant_to_size", {
  strategy <- FFW16_Strategy()
  pl <- FFW16_PlantPlus(strategy)
  sizes <- c(1, 5, 10, 12, strategy$hmat)
  env <- fixed_environment(1.0)
  res <- grow_plant_to_size(pl, sizes, "height", env, 10000)

  expect_that(res$state[, "height"], equals(sizes, tolerance=1e-6))

  sizes2 <- c(sizes, last(sizes) * 2)
  expect_that(res2 <- grow_plant_to_size(pl, sizes2, "height", env, 100),
              gives_warning("Time exceeded time_max"))
  expect_that(length(res2$time), equals(length(sizes2)))
  expect_that(last(res2$time), equals(NA_real_))
  expect_that(any(is.na(res2$time[-length(sizes2)])), is_false())

  expect_that(res3 <- grow_plant_to_size(pl, sizes2, "height", env,
                                         100, warn=FALSE),
              not(gives_warning()))
  expect_that(res3, equals(res2))

  expect_that(res4 <- grow_plant_to_size(pl, sizes2, "height", env,
                                         100, warn=FALSE, filter=TRUE),
              not(gives_warning()))

  ## Manually filter:
  cmp <- res2
  i <- !is.na(cmp$time)
  expect_that(res4$time, equals(cmp$time[i]))
  expect_that(res4$plant, equals(cmp$plant[i]))
  expect_that(res4$state, equals(cmp$state[i,]))
  expect_that(res4$trajectory, equals(cmp$trajectory))

  if (FALSE) {
    plot(height ~ time, as.data.frame(res$trajectory), type="l")
    points(res$time, res$state[, "height"], pch=19)

    plot(height ~ time, as.data.frame(res2$trajectory), type="l")
    points(res2$time, res2$state[, "height"], pch=19)
  }
})

test_that("grow_plant_to_time", {
  strategy <- FFW16_Strategy()
  pl <- FFW16_PlantPlus(strategy)
  env <- fixed_environment(1.0)
  times <- c(0, 10^(-4:3))
  res <- grow_plant_to_time(pl, times, env)
  expect_that(res$plant, is_a("list"))
  expect_that(length(res$plant), equals(length(times)))

  expect_that(res$state, is_a("matrix"))
  expect_that(colnames(res$state), equals(pl$ode_names))
  expect_that(nrow(res$state), equals(length(times)))

  expect_that(all(diff(res$state[, "height"]) > 0), is_true())

  expect_that(res$time, is_identical_to(times))
})
