## Support functions for growing a plant to a given size, etc.  Used
## in dfalster/TraitGrowthTrajectories

##' Grow a plant up to particular sizes.
##'
##' @title Grow plant to given size
##' @param plant A \code{Plant} object.
##' @param sizes A vector of sizes to grow the plant to (increasing in
##' size).
##' @param size_name The name of the size variable within
##' \code{Plant$rates} (e.g., height).
##' @param env An \code{Environment} object.
##' @param time_max Time to run the ODE out for -- only exists to
##' prevent an infinite loop (say, on an unreachable size).
##' @param warn Warn if requesting a plant that is too large?
##' @param filter Filter plants that are too large?
##' @return A list with elements \code{time} (the time that a given
##' size was reached), \code{state} (the \emph{ode state} at these
##' times, as a matrix) and \code{plant} a list of plants grown to the
##' appropriate size.  Note that if only a single size is given,
##' a list of length 1 is returned.
##' @export
grow_plant_to_size <- function(plant, sizes, size_name, env,
                               time_max=Inf, warn=TRUE, filter=FALSE) {
  obj <- grow_plant_bracket(plant, sizes, size_name, env, time_max, warn)

  polish <- function(i) {
    grow_plant_bisect(obj$runner, sizes[[i]], size_name,
                      obj$t0[[i]], obj$t1[[i]], obj$y0[i, ])
  }
  res <- lapply(seq_along(sizes), polish)

  state <- t(sapply(res, "[[", "state"))
  colnames(state) <- colnames(obj$state)

  ret <- list(time=vnapply(res, "[[", "time"),
              state=state,
              plant=lapply(res, "[[", "plant"),
              trajectory=cbind(time=obj$time, state=obj$state),
              env=env)
  if (filter) {
    i <- !vlapply(ret$plant, is.null)
    if (!all(i)) {
      ret$time  <- ret$time[i]
      ret$state <- ret$state[i, , drop=FALSE]
      ret$plant <- ret$plant[i]
    }
  }
  ret
}

##' @export
##' @rdname grow_plant_to_size
##' @param ... Additional parameters passed to
##' \code{grow_plant_to_size}.
##' @param heights Heights (when using \code{grow_plant_to_height})
grow_plant_to_height <- function(plant, heights, env, ...) {
  grow_plant_to_size(plant, heights, "height", env, ...)
}

##' Grow a plant up for particular time lengths
##'
##' @title Grow a plant
##' @param plant A \code{Plant} object
##' @param times A vector of times
##' @param env An \code{Environment} object
##' @export
grow_plant_to_time <- function(plant, times, env) {
  if (any(times < 0.0)) {
    stop("Times must be positive")
  }
  n <- length(times)
  if (n == 0L) {
    stop("At least one time must be given")
  }

  y0 <- plant$ode_state
  t0 <- 0.0
  i <- 1L
  t_next <- times[[i]]
  strategy_name <- plant$strategy_name

  pr1 <- PlantRunner(strategy_name, "Env")(plant, env)
  pr2 <- PlantRunner(strategy_name, "Env")(plant, env)

  runner <- OdeRunner(strategy_name)(pr1)
  runner_detail <- OdeRunner(strategy_name)(pr2)

  ## TODO: This could also be done by better configuring the
  ## underlying ODE runner, but this seems a reasonable way of getting
  ## things run for now.
  state <- matrix(NA, n, length(y0))
  colnames(state) <- plant$ode_names

  plant <- vector("list", n)
  while (i <= n) {
    runner$step()
    t1 <- runner$time
    y1 <- runner$state
    while (t_next < t1 && i <= n) {
      runner_detail$set_state(y0, t0)
      runner_detail$step_to(t_next)
      state[i, ] <- runner_detail$state
      plant[[i]] <- runner_detail$object$plant
      i <- i + 1L
      t_next <- times[i] # allows out-of-bounds extraction
    }
    t0 <- t1
    y0 <- y1
  }

  list(time=times, state=state, plant=plant, env=env)
}

## internal funciton to grab the internal state of the ode runner:
get_plant_internals_fun <- function (plant) {
  get(paste0(plant$strategy_name, '_oderunner_plant_internals'))
}

grow_plant_bracket <- function(plant, sizes, size_name, env,
                               time_max=Inf, warn=TRUE) {
  if (length(sizes) == 0L || is.unsorted(sizes)) {
    stop("sizes must be non-empty and sorted")
  }
  if (plant$state(size_name) > sizes[[1]]) {
    stop("Plant already bigger than smallest target size")
  }
  strategy_name <- plant$strategy_name

  # TODO: size index uses index from 0
  # can we clarify?
  size_index <- (which(plant$ode_names == size_name) - 1)

  runner <- OdeRunner(strategy_name)(PlantRunner(strategy_name, "Env")(plant, env))
  internals <- get_plant_internals_fun(runner$object$plant)
  i <- 1L
  n <- length(sizes)
  j <- rep_len(NA_integer_, n)
  state <- list(list(time=runner$time, state=runner$state))

  while (i <= n & runner$time < time_max) {
    ok <- tryCatch({
      runner$step()
      TRUE
    },
    error=function(e) {
      msg <- c("Stopping early as integration failed with error: ", e$message,
               sprintf("  %d larger sizes dropped", sum(is.na(j))))
      if (warn) {
        warning(paste(msg, collapse="\n"), immediate.=TRUE)
      }
      ## TODO: Consider making this an error, or making the test a bit better.
      if (runner$object$plant$ode_rates[[2]] < 1e-10) {
        warning("Integration may have failed for reasons other than mortality",
                immediate.=TRUE)
      }
      FALSE
    })
    if (!ok) {
      break
    }
    state <- c(state, list(list(time=runner$time, state=runner$state)))


    while (i <= n && internals(runner)$state(size_index) > sizes[[i]]) {
      j[[i]] <- length(state) - 1L
      i <- i + 1L
    }
    if (runner$time >= time_max) {
      if (warn) {
        warning(sprintf("Time exceeded time_max, %d larger sizes dropped",
                        sum(is.na(j))), immediate.=TRUE)
      }
      break
    }
  }

  t <- vnapply(state, "[[", "time")
  m <- t(sapply(state, "[[", "state"))
  k <- j + 1L
  colnames(m) <- runner$object$plant$ode_names
  list(t0=t[j],
       t1=t[k],
       y0=m[j, , drop=FALSE],
       y1=m[k, , drop=FALSE],
       time=t,
       state=m,
       index=j,
       runner=runner)
}

##' @noRd
##' @importFrom stats uniroot
grow_plant_bisect <- function(runner, size, size_name, t0, t1, y0) {

  
  # TODO: size index uses index from 0
  # can we clarify?
  size_index <- (which(runner$object$plant$ode_names == size_name) - 1)

  internals <- get_plant_internals_fun(runner$object$plant)
  f <- function(t1) {
    runner$set_state(y0, t0)
    runner$step_to(t1)
    internals(runner)$state(size_index) - size
  }

  if (is.na(t0) || is.na(t1) || any(is.na(y0))) {
    y0[] <- NA_real_
    list(time=NA_real_, state=y0, plant=NULL)
  } else {
    root <- uniroot(f, lower=t0, upper=t1)
    list(time=root$root, state=runner$state, plant=runner$object$plant)
  }
}

## These are waiting on RcppR6 #23 and plant #164

# ## This will get merged into RcppR6, so may change!
# plant_to_plant_plus <- function(x, ...) {
#   UseMethod("plant_to_plant_plus")
# }
# ##' @export
# `plant_to_plant_plus.Plant<FF16>` <- function(x, ...) {
#   FF16_plant_to_plant_plus(x, ...)
# }
# ##' @export
# `plant_to_plant_plus.Plant<FF16r>` <- function(x, ...) {
#   FF16r_plant_to_plant_plus(x, ...)
# }

#' Compute the whole plant light compensation point for a single
#' plant.
#' @title Whole plant light compensation point
#' @param p A \code{PlantPlus}, with strategy, height, etc set.
#' @param ... Additional arguments that are ignored
#' @export
#' @author Rich FitzJohn
lcp_whole_plant <- function(p, ...) {
  UseMethod("lcp_whole_plant")
}
##' @export
`lcp_whole_plant.Plant<FF16>` <- function(p, ...) {
  FF16_lcp_whole_plant(p, ...)
}
##' @export
`lcp_whole_plant.Plant<FF16r>` <- function(p, ...) {
  FF16r_lcp_whole_plant(p, ...)
}
##' @export
lcp_whole_plant.Plant <- function(p, ...) {
  lcp_whole_plant(p, ...)
}
