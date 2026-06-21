#' Turn `species` component of plant solver output into a tidy data object 
#'
#' @rdname tidy_patch
#' @importFrom rlang .data
tidy_species <- function(results) {

  n_spp <- length(results[[1]]$species)

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
  cumulative_names <- c("sum_rainfall","sum_infiltration","sum_drainage","sum_resource_depletion")
  
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


#' Interpolate data on size distributions for each species to specific timer points, using specified interpolation method
#'
#' @param tidy_species_data output of either `tidy_patch` or `tidy_species`
#' @param times times to interpolate to
#' @param method Method for interpolation. For more info see help on stats::spline
#'
#' @return Returns a tibble of similar format to input, but with all outputs interpolated to specified hieghts.
#' @export
#' @importFrom stats spline
#' @importFrom rlang .data
interpolate_to_times <- function(tidy_species_data, times, method="natural") {
  
  # helper function - predicts to new values with spline
  # only needed to ensure predictions for xout outside the ange of x are set to NA
  f <- function(x, y, xout) {
    y_pred <- stats::spline(x, y, xout=xout, method=method)$y
    y_pred[!dplyr::between(xout, min(x), max(x))] <- NA
    y_pred
  }
  
  tidy_species_data %>%
    tidyr::drop_na() %>%
    dplyr::group_by(.data$species, .data$node) %>%
    dplyr::reframe(
      dplyr::across(tidyselect::where(is.double), ~f(.data$time, .x, xout=times)),
    )
}


#' Interpolate data on size distributions for each species to specific heights at every time point
#'
#' @param tidy_species_data output of either `tidy_patch` or `tidy_species`
#' @param heights heights to interpolate to
#' @param method Method for interpolation. For more info see help on stats::spline
#' @param min_log_density Set minimum possible value of log_density
#' @return Returns a tibble of similar format to input, but with all outputs interpolated to specified hieghts.
#' @export
#' @importFrom stats spline
#' @importFrom rlang .data
interpolate_to_heights <- function(tidy_species_data, heights, method="natural", min_log_density = -100) {
  
  # helper function - predicts to new values with spline
  # only needed to ensure predictions for xout outside the ange of x are set to NA
  f <- function(x, y, xout, time) {

    y_pred <- stats::spline(x, y, xout=xout, method=method)$y
    y_pred[!dplyr::between(xout, min(x), max(x))] <- NA

    y_pred
  }

  if(!exists("step", tidy_species_data)){
    tidy_species_data <- tidy_species_data %>%
      tibble::add_column(step = NA)
  }
  
  tidy_species_data %>%
    tidyr::drop_na(-dplyr::any_of("step")) %>%
    
    # check for very negative values
    dplyr::mutate(
      log_density = ifelse(.data$log_density < min_log_density, min_log_density, .data$log_density)
    ) %>%
    dplyr::group_by(.data$species, .data$time, .data$step) %>%
    # remove any repated x values - these cuase warnings in the interpolation
    dplyr::filter(!duplicated(.data$height)) %>%

    dplyr::reframe(
      dplyr::across(tidyselect::where(is.double), ~ f(.data$height, .x, xout = heights, .data$time[1]))
    ) %>%
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
#' Integrations are performed using trapezium integration
#'
#' @param tidy_species_data output of either `tidy_patch` or `tidy_species`
#'
#' @return a tibble whose columns provide metrics on integrated totals for each variable for each species at each time

#' @export
#'
#' @importFrom rlang .data
integrate_over_size_distribution <- function(tidy_species_data) {
  tidy_species_data  %>%
    dplyr::select(-dplyr::any_of("node")) %>% stats::na.omit() %>%
    dplyr::filter(.data$step > 1) %>% 
    dplyr::group_by(.data$step, .data$time, .data$patch_density, .data$species) %>% 
    dplyr::reframe(
      density_integrated = -trapezium(.data$height, .data$density), 
      min_height = min(.data$height),
      dplyr::across(tidyselect::where(is.double) & !c(.data$density, .data$density_integrated, .data$min_height) , ~-trapezium(height, density * .x)) 
    ) %>% 
    dplyr::rename(density = .data$density_integrated)
}
