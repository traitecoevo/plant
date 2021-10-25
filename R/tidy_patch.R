

#' Turn `species` component of scm output into a tidy data object 
#'
#' @param data a list, the `species` component of scm output.
#'
#' @return a tibble describing columns of the output for each cohort in each species
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
      tidyr::pivot_longer(cols=starts_with("V"), names_to = "cohort") %>% 
      dplyr::pull(value)
  }
  
  data_species %>% dplyr::mutate(density = exp(log_density)) %>% dplyr::select(-log_density)
}



#' Turn `env` component of scm output into a tidy data object 
#'
#' @param data a list, the `env` component of scm output.
#'
#' @return a tibble describing the environment in a patch
tidy_env <- function(env) {
  dplyr::tibble(step = seq_len(length(env))) %>%
    dplyr::left_join(by = "step", 
    env %>% purrr::map_df(tidyr::as_tibble, .id= "step") %>% dplyr::mutate(step = as.integer(step))
    )
}





#' Turns scm output into a tidy data object 
#'
#' @param data output of run_scm_collect
#'
#' @return a list
#' @export
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


#' XXXXX
#'
#' @param data 
#'
#' @return
#' @export
#'
#' @examples
patch_species_total <- function(data) {
  data  %>%
    select(-cohort) %>% na.omit() %>% 
    filter(step > 1) %>% 
    group_by(step, time, patch_density, species) %>% 
    summarise(
      individuals = -plant:::trapezium(height, density),
      across(c(starts_with("area"), starts_with("mass")), ~ -plant:::trapezium(height, density*.x)),#, .names = "{.col}_tot"),
      .groups="drop"
    )
}
