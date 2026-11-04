#' Sample individuals from the size distribution
#'
#' Draws a random sample of individual plants from a stand simulated by
#' [run_scm()]: how many stand in a given ground area, and how big each one is.
#' Given several times, the same individuals are followed through all of them,
#' so each grows with its cohort and, once dead, stays dead.
#'
#' The SCM carries density per unit ground area and per unit birth date.
#' Individuals are drawn as a Poisson process over birth date, at the density at
#' the first requested time they could be seen, interpolated between nodes with
#' a shape-preserving cubic (PCHIP). Each then carries a fixed hazard budget, drawn from a
#' unit exponential, and is dead from the first requested time at which its
#' cumulative mortality since then exceeds that budget. In birth date a cohort's
#' density changes only through mortality, so at every requested time the number
#' of individuals is Poisson, with mean the area times the integral of the
#' density over birth date, and the sample is exact apart from the interpolation
#' between nodes. Every other state variable is interpolated to an individual's
#' birth date in the same way.
#'
#' The run must carry density in birth date,
#' `control(node_density_coordinate = "birth_date")`. Set the seed with
#' [set.seed()] to make a sample reproducible.
#'
#' @param results Output of `run_scm(p, ctrl, collect = TRUE)`, run with
#'   density in birth date.
#' @param area Ground area to sample (m^2).
#' @param times Times at which to report the sample. Defaults to the end of the
#'   run, a single snapshot.
#' @return A tibble with one row per individual alive at each requested time:
#'   `species`, `individual` (the same at every time), `birth_time`, `time`, and
#'   each node state variable interpolated to the individual (`height`,
#'   `mortality`, `fecundity`, ...).
#' @export
#' @examples
#' \dontrun{
#' p <- scm_base_parameters("FF16")
#' p$max_patch_lifetime <- 50
#' p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
#' results <- run_scm(p, ctrl = control(node_density_coordinate = "birth_date"),
#'                    collect = TRUE)
#' set.seed(1)
#' trees <- sample_size_distribution(results, area = 400, times = c(10, 20, 50))
#' }
sample_size_distribution <- function(results, area, times = NULL) {
  sp <- results$species
  if (!"log_density_state" %in% names(sp)) {
    stop("Sampling needs density in birth date; run with ",
         "control(node_density_coordinate = \"birth_date\")", call. = FALSE)
  }
  if (!(length(area) == 1 && is.finite(area) && area > 0)) {
    stop("'area' must be one positive number", call. = FALSE)
  }
  t_end <- max(sp$time)
  times <- if (is.null(times)) t_end else times
  if (!is.numeric(times) || length(times) == 0 || any(!is.finite(times))) {
    stop("'times' must be finite numbers", call. = FALSE)
  }
  times <- sort(unique(times))
  if (any(times < 0 | times > t_end)) {
    stop("'times' must lie within the run, from 0 to ", t_end, call. = FALSE)
  }
  columns <- state_columns(sp)

  out <- lapply(sort(unique(sp$species)), function(s) {
    knots <- birth_date_knots(sp[sp$species == s, ], times, columns)
    x <- sample_species(knots, times, area, columns)
    if (nrow(x) > 0) x$species <- s
    x
  })
  out <- dplyr::bind_rows(out)
  if (nrow(out) == 0) {
    ## Same columns whether or not anything was drawn.
    empty <- c(list(species = character(), individual = integer(),
                    birth_time = numeric(), time = numeric()),
               stats::setNames(rep(list(numeric()), length(columns)), columns))
    return(tibble::as_tibble(empty))
  }
  key <- paste(out$species, out$individual)
  out$individual <- match(key, unique(key))
  tibble::as_tibble(out[c("species", "individual", "birth_time", "time", columns)])
}

## Number per unit area at each time, integrated over birth date with the same
## interpolant the sampler draws from; a matrix of species by time.
birth_date_totals <- function(results, times) {
  sp <- results$species
  spp <- sort(unique(sp$species))
  out <- vapply(spp, function(s) {
    knots <- birth_date_knots(sp[sp$species == s, ], times, character(0))
    vapply(seq_along(times), function(j) {
      k <- knots[[j]]
      if (length(k$tau) < 2) return(0)
      sum(spline_segments(spline_density(k), k$tau))
    }, 0)
  }, numeric(length(times)))
  matrix(out, nrow = length(spp), byrow = TRUE, dimnames = list(spp, times))
}

## The node state variables to carry to individuals: everything numeric that is
## a property of the plant rather than of the density or the patch.
state_columns <- function(sp) {
  setdiff(names(sp)[vapply(sp, is.double, TRUE)],
          c("time", "patch_density", "density", "log_density",
            "log_density_state"))
}

## For each requested time, the nodes as knots in birth date: their birth dates
## `tau`, the density per unit birth date `n`, cumulative mortality `m` and the
## state columns, with the boundary (the cohort being born at that moment) as
## the last knot.
##
## Node k first appears as the boundary at its own birth date, so its first row
## gives that date. Between collected steps a node's state is interpolated in
## time with pchip() (height, mortality and heartwood only rise), and
## the boundary's from the boundary rows of the steps either side.
birth_date_knots <- function(d, times, columns) {
  d <- d[order(d$node, d$time), ]
  ## A node with no survivors carries an enormous or infinite hazard, which a
  ## spline cannot pass through; exp(-750) is zero either way.
  d$mortality <- pmin(d$mortality, 750)
  d$log_density_state <- pmax(d$log_density_state, -750)
  cols <- unique(c("log_density_state", "mortality", columns))
  by_node <- lapply(split(d, d$node), function(x) x[!duplicated(x$time), ])
  birth <- vapply(by_node, function(x) x$time[1], 0)

  at_times <- function(time, value) {
    out <- rep(NA_real_, length(times))
    ok <- times >= min(time) & times <= max(time)
    out[ok] <- if (length(time) == 1) value else pchip(time, value)(times[ok])
    out
  }
  ## times x nodes, one matrix per column
  node_vals <- lapply(stats::setNames(cols, cols), function(col) {
    matrix(vapply(by_node, function(x) at_times(x$time, x[[col]]),
                  numeric(length(times))),
           nrow = length(times))
  })
  first_rows <- do.call(rbind, lapply(by_node, function(x) x[1, ]))
  boundary <- lapply(stats::setNames(cols, cols), function(col) {
    if (nrow(first_rows) == 1) rep(first_rows[[col]], length(times)) else
      stats::approx(first_rows$time, first_rows[[col]], times, rule = 2)$y
  })

  lapply(seq_along(times), function(j) {
    born <- which(birth < times[j] & !is.na(node_vals$mortality[j, ]))
    k <- list(tau = c(unname(birth[born]), times[j]))
    for (col in cols) {
      k[[col]] <- c(node_vals[[col]][j, born], boundary[[col]][j])
    }
    k$n <- exp(k$log_density_state)
    k
  })
}

spline_density <- function(k) {
  pchip(k$tau, k$n)
}

## Piecewise cubic Hermite interpolation that keeps each segment within its two
## end values (Fritsch & Carlson 1980, with the slopes of Fritsch & Butland 1984,
## as MATLAB's pchip). stats::splinefun(method = "monoH.FC") promises this only
## for monotone data: across a cohort that failed to establish, mortality jumps
## to the 750 sentinel, and that spline then overshoots to negative values
## between ordinary knots. Returns f(x, deriv = 0) like splinefun().
pchip <- function(x, y) {
  n <- length(x)
  h <- diff(x)
  delta <- diff(y) / h
  d <- numeric(n)
  if (n == 2) {
    d[] <- delta
  } else {
    a <- delta[-(n - 1)]
    b <- delta[-1]
    w1 <- 2 * h[-1] + h[-(n - 1)]
    w2 <- h[-1] + 2 * h[-(n - 1)]
    d[2:(n - 1)] <- ifelse(a * b > 0, (w1 + w2) / (w1 / a + w2 / b), 0)
    end_slope <- function(h1, h2, d1, d2) {
      s <- ((2 * h1 + h2) * d1 - h1 * d2) / (h1 + h2)
      if (sign(s) != sign(d1)) 0
      else if (sign(d1) != sign(d2) && abs(s) > abs(3 * d1)) 3 * d1
      else s
    }
    d[1] <- end_slope(h[1], h[2], delta[1], delta[2])
    d[n] <- end_slope(h[n - 1], h[n - 2], delta[n - 1], delta[n - 2])
  }
  function(xout, deriv = 0) {
    i <- findInterval(xout, x, rightmost.closed = TRUE, all.inside = TRUE)
    s <- xout - x[i]
    c2 <- (3 * delta[i] - 2 * d[i] - d[i + 1]) / h[i]
    c3 <- (d[i] + d[i + 1] - 2 * delta[i]) / h[i]^2
    if (deriv == 0) y[i] + s * (d[i] + s * (c2 + s * c3))
    else d[i] + s * (2 * c2 + 3 * s * c3)
  }
}

## Integral of the interpolant over each segment between breakpoints `x` (which
## must include every knot inside them): exact for a cubic Hermite piece, from
## its end values and slopes.
spline_segments <- function(f, x) {
  y <- f(x)
  dy <- f(x, deriv = 1)
  h <- diff(x)
  pmax(h * (y[-length(y)] + y[-1]) / 2 + h^2 * (dy[-length(dy)] - dy[-1]) / 12, 0)
}

## Draw `n` birth dates from the spline density on the segments between
## breakpoints `x`, by segment and then by rejection: pchip() keeps each segment
## monotone, so its larger end value bounds it.
draw_birth_dates <- function(f, x, n) {
  seg <- sample.int(length(x) - 1, n, replace = TRUE,
                    prob = spline_segments(f, x))
  lo <- x[seg]
  hi <- x[seg + 1]
  top <- pmax(f(lo), f(hi))
  tau <- numeric(n)
  todo <- seq_len(n)
  while (length(todo) > 0) {
    try <- stats::runif(length(todo), lo[todo], hi[todo])
    keep <- stats::runif(length(todo)) * top[todo] <= f(try)
    tau[todo[keep]] <- try[keep]
    todo <- todo[!keep]
  }
  tau
}

## Individuals of one species, followed through the requested times.
sample_species <- function(knots, times, area, columns) {
  ## Birth dates in (times[j-1], times[j]] are first seen at times[j].
  born <- lapply(seq_along(times), function(j) {
    k <- knots[[j]]
    if (length(k$tau) < 2) return(numeric(0))
    ## Never below the first birth date: the interpolant would extrapolate a
    ## density where no cohort exists (a resumed run, or a late first birth).
    lower <- if (j == 1) min(k$tau) else max(times[j - 1], min(k$tau))
    x <- sort(unique(c(lower, k$tau[k$tau > lower])))
    if (length(x) < 2) return(numeric(0))
    f <- spline_density(k)
    n <- stats::rpois(1, area * sum(spline_segments(f, x)))
    draw_birth_dates(f, x, n)
  })
  first <- rep(seq_along(times), lengths(born))
  tau <- unlist(born)
  budget <- stats::rexp(length(tau))
  m_first <- rep(NA_real_, length(tau))
  alive <- rep(TRUE, length(tau))

  rows <- vector("list", length(times))
  for (j in seq_along(times)) {
    i <- which(first <= j & alive)
    if (length(i) == 0) next
    k <- knots[[j]]
    at <- function(col) {
      if (length(k$tau) < 2) rep(k[[col]], length(i)) else
        pchip(k$tau, k[[col]])(tau[i])
    }
    m <- at("mortality")
    new <- first[i] == j
    m_first[i[new]] <- m[new]
    dead <- m - m_first[i] > budget[i]
    alive[i[dead]] <- FALSE
    i <- i[!dead]
    if (length(i) == 0) next
    vals <- lapply(stats::setNames(columns, columns), at)
    rows[[j]] <- data.frame(individual = i, birth_time = tau[i], time = times[j],
                            vals, check.names = FALSE)
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) data.frame(individual = integer(), birth_time = numeric(),
                               time = numeric()) else out
}
