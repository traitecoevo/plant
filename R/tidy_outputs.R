#' Turn `species` component of plant solver output into a tidy data object 
#'
#' @rdname tidy_patch
#' @importFrom rlang .data
tidy_species <- function(results) {

  n_spp <- length(results[[1]]$species)
  warn_top_layer(results)

   get_species_sdd <- function(i) {
     purrr::imap_dfr(results, ~ .x$species[[i]] |>
       t() |>
       dplyr::as_tibble() |>
       dplyr::mutate(step = .y, node = seq_len(dplyr::n()), species = i))
   }
  
  purrr::map_dfr(seq_len(n_spp), get_species_sdd) |>
    dplyr::mutate(
      density = exp(.data$log_density),
      species = as.character(.data$species)
    ) 
}


## The absorbing top layer (control$absorb_epsilon) is not a node, so the rows
## tidy_species() returns, and every integral over them, leave out what it
## carries. Said rather than refused: this runs at the end of a collected run.
warn_top_layer <- function(results) {
  top <- unlist(lapply(results, function(r)
    lapply(r$species, function(s) attr(s, "top_density"))))
  if (length(top) > 0 && any(top > 0)) {
    warning("The absorbing top layer carries density (up to ",
            signif(max(top), 3), ") that these node rows leave out, so ",
            "integrals over them undercount the stand. Read the layer from ",
            "species$top_layer or tidy_snapshots().", call. = FALSE)
  }
}

#' Turn `env` component of solver output into a tidy data object
#'
#' @rdname tidy_patch
#' @importFrom rlang .data
tidy_env <- function(results) {
  env <- lapply(results, "[[", "env")

  # get list of variables
  env_variables = names(env[[1]])
  
  # extract over each variable, concatenating by time, 
  # then join all variables by step, and force unnamed vectors 
  # to have variable names using regex magic
  env_long <- env_variables %>%
         purrr::map(
          function(v) purrr::map_dfr(env, 
            ~purrr::pluck(.x, v) %>% 
              data.frame, .id = "step") %>%
              dplyr::as_tibble() %>%
              dplyr::mutate(dplyr::across(dplyr::any_of("step"), as.integer)) %>%
              dplyr::rename_with(~ gsub("\\.", v, .x)
                                 
                                 
          )
        )

  names(env_long) <- env_variables
  if(any(env_variables == "soil_moist_cumulative_flux")){
  ## Positional: must stay in lockstep with the trailing state slots written in
  ## tf24_environment.h's compute_rates()/add_water_pulse().
  cumulative_names <- c("sum_rainfall","sum_infiltration","sum_drainage",
                        "sum_resource_depletion","sum_pulse_runoff")
  
  env_long$soil_moist_cumulative_flux <- env_long$soil_moist_cumulative_flux %>%
    dplyr::mutate(cumulative_variables = rep(cumulative_names, times = dplyr::n()/length(cumulative_names))) %>%
    tidyr::pivot_wider(names_from = "cumulative_variables", values_from = "soil_moist_cumulative_flux")
  }
  return(env_long)
}


#' Turns output of plant solver into a tidy data object 
#'
#' @param results output of run_scm(..., collect = TRUE)
#'
#' @return a list, containing outputs of plant solver in tidy format
#' @keywords internal
#' @importFrom rlang .data
tidy_patch <- function(results) {
  time <- sapply(results, "[[", "time")
  patch_density <- sapply(results, "[[", "patch_density")

  out <- list()

  out[["steps"]] <-
    dplyr::tibble(
      step = seq_len(length(time)),
      time = time,
      patch_density = patch_density
    )

  out[["n_spp"]] <- length(results[[1]]$species)

  out[["species"]] <-
    results |>
    tidy_species() |>
    dplyr::left_join(by = "step", out[["steps"]]) |>
    dplyr::select(dplyr::all_of(c("species", "time", "step", "patch_density", "node", "density", "log_density")), dplyr::everything())

  out[["env"]] <- 
    results |>
    tidy_env() |>
    purrr::map(dplyr::left_join, out[["steps"]], by = "step") |>
    purrr::map(~.x |> dplyr::select(dplyr::all_of(c("time", "step", "patch_density")), dplyr::everything()))
  
  out
}


#' Interpolate data on size distributions for each species to specific time points
#'
#' Each node's state is interpolated through time with a shape-preserving cubic
#' (PCHIP) by default, which stays between the values at neighbouring steps: a
#' height or a cumulative mortality never falls, and a density never goes
#' negative. A natural spline overshoots wherever a node's trajectory bends
#' sharply, as when a cohort fails to establish.
#'
#' @param tidy_species_data output of either `tidy_patch` or `tidy_species`
#' @param times times to interpolate to
#' @param method `"pchip"` (the default), or a method for [stats::spline()],
#'   such as `"natural"`, which this function used before.
#'
#' @return Returns a tibble of similar format to input, but with all outputs interpolated to specified times.
#' @export
#' @importFrom stats spline
#' @importFrom rlang .data
interpolate_to_times <- function(tidy_species_data, times, method = "pchip") {
  f <- interpolator(method)
  tidy_species_data %>%
    tidyr::drop_na() %>%
    dplyr::group_by(.data$species, .data$node) %>%
    dplyr::reframe(
      dplyr::across(tidyselect::where(is.double), ~f(.data$time, .x, xout=times)),
    )
}

## Interpolate y(x) to xout, NA outside the range of x. "pchip" keeps every
## segment between its end values; anything else is a stats::spline() method.
interpolator <- function(method) {
  function(x, y, xout) {
    keep <- !duplicated(x)
    x <- x[keep]
    y <- y[keep]
    o <- order(x)
    x <- x[o]
    y <- y[o]
    out <- rep(NA_real_, length(xout))
    ok <- !is.na(xout) & xout >= x[1] & xout <= x[length(x)]
    if (!any(ok)) return(out)
    out[ok] <- if (length(x) == 1) y else if (method == "pchip")
      pchip(x, y)(xout[ok]) else stats::spline(x, y, xout = xout[ok], method = method)$y
    out
  }
}


#' Interpolate data on size distributions for each species to specific heights at every time point
#'
#' At each time, the state is interpolated across nodes in height with a
#' shape-preserving cubic (PCHIP) by default, and density through its logarithm,
#' so it stays between the values at neighbouring nodes. Density can fall several
#' hundredfold between neighbouring nodes near the seedling end, and a natural
#' spline through its logarithm then overshoots by orders of magnitude.
#'
#' @param tidy_species_data output of either `tidy_patch` or `tidy_species`
#' @param heights heights to interpolate to
#' @param method `"pchip"` (the default), or a method for [stats::spline()],
#'   such as `"natural"`, which this function used before.
#' @param min_log_density Set minimum possible value of log_density
#' @return Returns a tibble of similar format to input, but with all outputs interpolated to specified hieghts.
#' @export
#' @importFrom stats spline
#' @importFrom rlang .data
interpolate_to_heights <- function(tidy_species_data, heights, method = "pchip", min_log_density = -100) {

  f <- interpolator(method)

  if(!exists("step", tidy_species_data)){
    tidy_species_data <- tidy_species_data %>%
      tibble::add_column(step = NA)
  }

  prepared <- tidy_species_data %>%
    tidyr::drop_na(-dplyr::any_of("step")) %>%
    # check for very negative values
    dplyr::mutate(
      log_density = ifelse(.data$log_density < min_log_density, min_log_density, .data$log_density)
    ) %>%
    dplyr::group_by(.data$species, .data$time, .data$step) %>%
    # remove any repeated x values - these cause warnings in the interpolation
    dplyr::filter(!duplicated(.data$height))

  # interpolate each state onto the requested heights; grid points outside the
  # observed range interpolate to NA and are dropped
  interpolated <- prepared %>%
    dplyr::reframe(
      dplyr::across(tidyselect::where(is.double), ~ f(.data$height, .x, xout = heights))
    ) %>%
    tidyr::drop_na(dplyr::any_of("height"))

  # On a coarse grid the tallest node usually sits above the highest in-range
  # grid point, so it would be dropped and the largest size bracket lost.
  # Append the actual largest individual in each time step instead, rather than
  # misrepresenting its density at the highest interpolation point (#352).
  largest <- prepared %>%
    dplyr::filter(.data$height == max(.data$height)) %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::any_of(names(interpolated)))

  dplyr::bind_rows(interpolated, largest) %>%
    dplyr::arrange(.data$species, .data$time, .data$step, .data$height) %>%
    dplyr::mutate(density = exp(.data$log_density))
}


#' Turn `results` of plant solver, when solving individuals into a tidy data object
#'
#' @param results plant solver output.
#'
#' @return a tibble whose columns provide metrics on each individual over time
#'
#' @export
tidy_individual <- function(results) {
  out <- dplyr::tibble(
    step = seq_len(length(results$time)),
    time = results$time
  ) %>%
    dplyr::bind_cols(
      height = results$state[, 1],
      mortality = results$state[, 2],
      fecundity = results$state[, 3],
      area_heartwood = results$state[, 4],
      mass_heartwood = results$state[, 5]
    )

  out
}


#' Integrate over the size distribution for each species at each time point, to give totals of each variable
#'
#' A run carried in birth date (`control(node_density_coordinate =
#' "birth_date")`, the default for the TF24 models) is integrated over birth
#' date: the density per unit birth date, and its product with each variable,
#' are interpolated between nodes with a shape-preserving cubic and integrated
#' exactly. Across a cohort that has just formed, density per unit height can
#' fall several hundredfold, and a trapezium in height overcounts there; in
#' birth date the integrand is smooth. A run carried in height is integrated over
#' height with the trapezium rule.
#'
#' @param tidy_species_data output of either `tidy_patch` or `tidy_species`.
#'   Integrating over birth date needs every step of the run, because a node's
#'   birth date is the time it first appears.
#' @param lower,upper heights bounding the integral; the defaults integrate over
#'   the whole size distribution. Over birth date, a bound cuts where the
#'   interpolated height of the cohorts crosses it. Over height, a bound that
#'   falls between two nodes takes the integrand by linear interpolation, as the
#'   trapezium rule already assumes. Either way, totals over adjacent ranges sum
#'   to the total over their union. For a finer height grid, interpolate with
#'   `interpolate_to_heights()` first.
#' @param coordinate `"auto"` (the default) integrates over birth date when the
#'   run carries density in birth date, and over height otherwise; `"birth_date"`
#'   and `"height"` choose one. `"height"` gives the trapezium totals this
#'   function returned before it integrated over birth date.
#'
#' @return a tibble whose columns provide metrics on integrated totals for each variable for each species at each time.
#'   `min_height` is the lowest height in the integrated range, and is `NA`
#'   when the range holds no individuals (the totals are then zero).

#' @export
#'
#' @importFrom rlang .data
integrate_over_size_distribution <- function(tidy_species_data,
                                             lower = -Inf, upper = Inf,
                                             coordinate = c("auto", "birth_date", "height")) {
  if (!(length(lower) == 1 && length(upper) == 1 && lower < upper)) {
    stop("'lower' and 'upper' must be single heights with lower < upper")
  }
  coordinate <- match.arg(coordinate)
  in_birth_date <- "log_density_state" %in% names(tidy_species_data)
  if (coordinate == "auto") {
    coordinate <- if (in_birth_date) "birth_date" else "height"
  }
  if (coordinate == "birth_date") {
    if (!in_birth_date) {
      stop("Integrating over birth date needs a run carried in birth date, ",
           "control(node_density_coordinate = \"birth_date\"); ",
           "or use coordinate = \"height\"", call. = FALSE)
    }
    return(integrate_over_birth_date(tidy_species_data, lower, upper))
  }
  tidy_species_data  %>%
    dplyr::select(-dplyr::any_of("node")) %>% stats::na.omit() %>%
    dplyr::filter(.data$step > 1) %>% 
    dplyr::group_by(.data$step, .data$time, .data$patch_density, .data$species) %>%
    # This one really does integrate over height, so a crossed grid is a live
    # hazard: reserve-gated growth lets a younger cohort overtake an older one,
    # and neighbouring trapezia then cancel instead of accumulating. Ordering the
    # grid first is what makes the rule below hold however the rows arrive.
    dplyr::arrange(.data$height, .by_group = TRUE) %>%
    dplyr::reframe(
      density_integrated = trapezium_between(.data$height, .data$density, lower, upper),
      min_height = if (lower < max(.data$height) && upper > min(.data$height)) max(lower, min(.data$height)) else NA_real_,
      dplyr::across(tidyselect::where(is.double) & !c("density", "density_integrated", "min_height") , ~trapezium_between(height, density * .x, lower, upper))
    ) %>%
    dplyr::rename(density = "density_integrated")
}

## The same totals, over birth date. Nodes are ordered by birth date, which
## never crosses, so unlike the height path the grid needs no reordering to be
## valid; a younger cohort overtaking an older one only makes the height bound
## cut a birth-date segment more than once, which birth_date_between() allows.
integrate_over_birth_date <- function(tidy_species_data, lower, upper) {
  tidy_species_data %>%
    add_birth_dates() %>%
    dplyr::select(-dplyr::any_of("node")) %>%
    dplyr::filter(.data$step > 1, !is.na(.data$log_density_state),
                  !is.na(.data$height)) %>%
    dplyr::group_by(.data$step, .data$time, .data$patch_density, .data$species) %>%
    dplyr::arrange(.data$birth_time, .by_group = TRUE) %>%
    dplyr::reframe(
      density_integrated = birth_date_between(.data$birth_time, .data$height,
                                              exp(.data$log_density_state),
                                              lower, upper),
      min_height = if (lower < max(.data$height) && upper > min(.data$height)) max(lower, min(.data$height)) else NA_real_,
      dplyr::across(tidyselect::where(is.double) &
                      !c("density", "density_integrated", "min_height", "birth_time"),
                    ~birth_date_between(birth_time, height,
                                        exp(log_density_state) * .x, lower, upper))
    ) %>%
    dplyr::rename(density = "density_integrated")
}

## A node's birth date is the time it first appears: it enters the collected
## output as the boundary at the moment it is born. Nodes seeded together at the
## start of a run (an initial state) all first appear at once, and then their
## birth dates cannot be told apart.
add_birth_dates <- function(tidy_species_data) {
  d <- tidy_species_data %>%
    dplyr::group_by(.data$species, .data$node) %>%
    dplyr::mutate(birth_time = min(.data$time)) %>%
    dplyr::ungroup()
  births <- unique(d[c("species", "node", "birth_time")])
  if (anyDuplicated(births[c("species", "birth_time")]) > 0) {
    stop("Nodes share a first appearance, so their birth dates cannot be ",
         "recovered: nodes seeded by an initial state, or a table filtered to ",
         "fewer steps than the run collected? Use coordinate = \"height\"",
         call. = FALSE)
  }
  d
}

## Integral over birth date of y, interpolated with pchip(), over the cohorts
## whose interpolated height lies in [lower, upper]. Height is monotone on each
## segment of its pchip(), so a bound crosses a segment at most once.
birth_date_between <- function(tau, h, y, lower, upper) {
  if (length(tau) < 2) return(0)
  ## A variable missing at a node has no total (log_density, the density per
  ## unit height, is missing wherever height does not vary with birth date).
  if (any(!is.finite(y))) return(NA_real_)
  f <- pchip(tau, y)
  if (lower <= min(h) && upper >= max(h)) {
    return(sum(spline_segments(f, tau)))
  }
  g <- pchip(tau, h)
  x <- tau
  for (b in c(lower, upper)) {
    if (!is.finite(b)) next
    for (i in which((h[-length(h)] - b) * (h[-1] - b) < 0)) {
      x <- c(x, stats::uniroot(function(t) g(t) - b, c(tau[i], tau[i + 1]),
                               tol = 1e-12 * (tau[i + 1] - tau[i]))$root)
    }
  }
  x <- sort(unique(x))
  mid <- g((x[-1] + x[-length(x)]) / 2)
  ## Half-open, [lower, upper), so that cohorts tied at a bound fall in exactly
  ## one of two adjacent ranges; closed above only at the tallest cohort.
  inside <- mid >= lower & (mid < upper | upper >= max(h))
  sum(spline_segments(f, x)[inside])
}

# Trapezium integral of y over sorted x, restricted to [lower, upper]. Each bound
# takes y by linear interpolation, the integrand the rule already assumes, so
# integrals over adjacent ranges add exactly. Bounds spanning x reduce to the
# plain rule.
trapezium_between <- function(x, y, lower, upper) {
  n <- length(x)
  if (lower <= x[1] && upper >= x[n]) {
    return(trapezium(x, y))
  }
  lower <- max(lower, x[1])
  upper <- min(upper, x[n])
  if (lower >= upper) {
    return(0)
  }
  at <- function(h) stats::approx(x, y, xout = h, ties = "ordered")$y
  inside <- x > lower & x < upper
  trapezium(c(lower, x[inside], upper), c(at(lower), y[inside], at(upper)))
}
