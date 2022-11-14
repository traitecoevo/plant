#' Turn `species` component of plant solver output into a tidy data object 
#'
#' @param data a list, the `species` component of plant solver output.
#'
#' @return a tibble whose columns provide metrics on each breakpoint in species size distribution
#' @importFrom rlang .data
tidy_species <- function(data) {
  
  # get dimensions of data = number of steps * number of nodes
  dimensions <- dim(data[1,,] )
  
  # establish data structure for results    
  data_species <- 
    tidyr::expand_grid(
      step = seq_len(dimensions[1]), 
      node = seq_len(dimensions[2])
    )
  
  # retrieve bnames of all tracked variables
  vars <- data[,1,1] %>% names()
  
  # bind each onto main data frame
  for(v in vars) {
    data_species[[v]] <- 
      data[v, , ] %>% 
      as.data.frame %>% tidyr::as_tibble() %>%
      tidyr::pivot_longer(cols=dplyr::starts_with("V"), names_to = "node") %>%
      dplyr::pull(.data$value)
  }
  
  data_species %>% dplyr::mutate(density = exp(.data$log_density))
}


#' Turn `env` component of solver output into a tidy data object 
#'
#' @param env a list, the `env` component of solver output.
#'
#' @return a tibble describing the environment in a patch
#' @importFrom rlang .data
tidy_env <- function(env) {
  # get list of variables
  env_variables = names(env[[1]])
  
  # extract over each variable, concatenating by time, 
  # then join all variables by step, and force unnamed vectors 
  # to have variable names using regex magic
  env_long <- env_variables %>%
    purrr::map(., function(v) purrr::map_dfr(env, ~ purrr::pluck(., v) %>% 
                                        data.frame, .id = "step") %>%
                 dplyr::mutate(dplyr::across(step, as.integer)) %>%
                 dplyr::rename_with(~ gsub("\\.", v, .)) %>%
                 tibble::as_tibble()
                 
                 ) 
  
  names(env_long) <- env_variables
  return(env_long)
}


#' Turns output of plant solver into a tidy data object 
#'
#' @param results output of run_scm_collect
#'
#' @return a list, containing outputs of plant solver in tidy format
#' @export
#' @importFrom rlang .data
tidy_patch <- function(results) {

  out <- results
  
  data <- dplyr::tibble(
    step = seq_len(length(results$time)),
    time = results$time, 
    patch_density = results$patch_density
    )
  
  out[["species"]] <- 
    dplyr::left_join(by = "step", data,
      purrr::map_df(results$species, tidy_species, .id="species")
    )
  
  out[["env"]] <- 
    tidy_env(results$env) %>%
    purrr::map(., dplyr::left_join, data, by = "step")
  
  out[["n_spp"]] <- length(results$species)
  
  out[["patch_density"]] <- NULL
  
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
    dplyr::summarise(
      dplyr::across(where(is.double), ~f(.data$time, .x, xout=times)),
      .groups = "keep") %>%
    dplyr::ungroup()
}


#' Interpolate data on size distributions for each species to specific heights at every time point
#'
#' @param tidy_species_data output of either `tidy_patch` or `tidy_species`
#' @param heights heights to interpolate to
#' @param method Method for interpolation. For more info see help on stats::spline
#'
#' @return Returns a tibble of similar format to input, but with all outputs interpolated to specified hieghts.
#' @export
#' @importFrom stats spline
#' @importFrom rlang .data
interpolate_to_heights <- function(tidy_species_data, heights, method="natural") {
  
  # helper function - predicts to new values with spline
  # only needed to ensure predictions for xout outside the ange of x are set to NA
  f <- function(x, y, xout) {
    y_pred <- stats::spline(x, y, xout=xout, method=method)$y
    y_pred[!dplyr::between(xout, min(x), max(x))] <- NA
    y_pred
  }

  if(!exists("step", tidy_species_data)){
    tidy_species_data <- tidy_species_data %>%
      tibble::add_column(step = NA)
  }

  tidy_species_data %>%
    tidyr::drop_na(-.data$step) %>%
    dplyr::group_by(.data$species, .data$time, .data$step) %>%
    dplyr::summarise(
      dplyr::across(where(is.double), ~f(.data$height, .x, xout=heights)),
      .groups = "keep") %>%
    dplyr::mutate(density = exp(.data$log_density)) %>%
    dplyr::ungroup()
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
    dplyr::select(-.data$node) %>% stats::na.omit() %>%
    dplyr::filter(.data$step > 1) %>% 
    dplyr::group_by(.data$step, .data$time, .data$patch_density, .data$species) %>% 
    dplyr::summarise(
      density_integrated = -trapezium(.data$height, .data$density), 
      min_height = min(.data$height),
      dplyr::across(where(is.double) & !c(.data$density, .data$density_integrated, .data$min_height) , ~-trapezium(height, density * .x)), 
      .groups = "drop"
    ) %>% 
    dplyr::rename(density = density_integrated)
}
