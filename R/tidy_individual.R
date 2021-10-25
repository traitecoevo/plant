#' Title
#'
#' @param results 
#'
#' @return
#' @export
#'
#' @examples
tidy_individual <- function(results) {
  
  out <- dplyr::tibble(
    step = seq_len(length(results$time)),
    time = results$time) %>%
    bind_rows(height = results$state[,1],
              mortality = results$state[,2],
              fecundity = results$state[,3],
              area_heartwood = results$state[,4],
              mass_heartwood = results$state[,5])
  
  data <- tibble(
    step = seq_len(length(results$time)),
    time = results$time, 
    patch_density = results$patch_density
  )
  
  out[["species"]] <- 
    left_join(by = "step", data,
              purrr::map_df(results$species, tidy_species, .id="species")
    )
  
  out[["env"]] <- 
    left_join(by = "step", data,
              tidy_env(results$env)
    )
  
  out[["n_spp"]] <- length(results$species)
  
  out[["patch_density"]] <- NULL
  
  out
}