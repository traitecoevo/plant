

#' Turn `species` component of scm output into a tidy data object 
#'
#' @param data a list, the `species` component of scm output.
#'
#' @return a tibble describing columns of the output for each cohort in each species
#' @importFrom rlang .data
tidy_species <- function(data) {
  
  # get dimensions of data = number of steps * number of cohorts
  dimensions <- dim(data[1,,] )
  
  # establish data structure for results    
  data_species <- 
    tidyr::expand_grid(
      step = seq_len(dimensions[1]), 
      cohort = seq_len(dimensions[2])
    )
  
  # retrieve bnames of all tracked variables
  vars <- data[,1,1] %>% names()
  
  # bind each onto main data frame
  for(v in vars) {
    data_species[[v]] <- 
      data[v, , ] %>% 
      as.data.frame %>% tidyr::as_tibble() %>%
      tidyr::pivot_longer(cols=dplyr::starts_with("V"), names_to = "cohort") %>% 
      dplyr::pull(.data$value)
  }
  
  data_species %>% dplyr::mutate(density = exp(.data$log_density))
}



#' Turn `env` component of scm output into a tidy data object 
#'
#' @param data a list, the `env` component of scm output.
#'
#' @return a tibble describing the environment in a patch
#' @importFrom rlang .data
tidy_env <- function(env) {
  dplyr::tibble(step = seq_len(length(env))) %>%
    dplyr::left_join(by = "step", 
    env %>% purrr::map_df(tidyr::as_tibble, .id= "step") %>% dplyr::mutate(step = as.integer(.data$step))
    )
}





#' Turns scm output into a tidy data object 
#'
#' @param data output of run_scm_collect
#'
#' @return a list
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
    dplyr::left_join(by = "step", data,
      tidy_env(results$env)
    )
  
  out[["n_spp"]] <- length(results$species)
  
  out[["patch_density"]] <- NULL
  
  out
}

#' Interpolate tidyspecies data to specific time points
#'
#' @param tidy_species_data 
#' @param times times to interpolate to
#' @param method Method for interpolation. For more info see help on stats::spline
#'
#' @return ??
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
    dplyr::group_by(.data$species, .data$cohort) %>%
    dplyr::summarise(
      dplyr::across(where(is.double), ~f(.data$time, .x, xout=times)),
      .groups = "keep") %>%
    dplyr::mutate(time=.data$times) %>%
    dplyr::ungroup()
}


#' Interpolate tidyspecies data to specific heights at every time point
#'
#' @param tidy_species_data ???
#' @param heights heights to interpolate to
#' @param method Method for interpolation. For more info see help on stats::spline
#'
#' @return ??
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
    dplyr::mutate(height=.data$heights,
                  density = exp(.data$log_density)) %>%
    dplyr::ungroup()
}



#' XXXXX
#'
#' @param data ???
#'
#' @return ??
#' @export
#'
#' @importFrom rlang .data
patch_species_total <- function(data) {
  data  %>%
    dplyr::select(-.data$cohort) %>% stats::na.omit() %>% 
    dplyr::filter(.data$step > 1) %>% 
    dplyr::group_by(.data$step, .data$time, .data$patch_density, .data$species) %>% 
    dplyr::summarise(
      individuals = -trapezium(.data$height, .data$density),
      min_h = min(.data$height),
      dplyr::across(c(dplyr::starts_with("area"), dplyr::starts_with("mass")), ~ -trapezium(height, density*.x)),#, .names = "{.col}_tot"),
      .groups="drop"
    )
}
