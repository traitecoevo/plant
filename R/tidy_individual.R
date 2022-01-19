#' Title
#'
#' @param results ???
#'
#' @return ??
#' @export
tidy_individual <- function(results) {
  
  out <- dplyr::tibble(
    step = seq_len(length(results$time)),
    time = results$time) %>%
    dplyr::bind_cols(height = results$state[,1],
              mortality = results$state[,2],
              fecundity = results$state[,3],
              area_heartwood = results$state[,4],
              mass_heartwood = results$state[,5])

  out
}