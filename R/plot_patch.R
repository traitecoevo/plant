

#' Title
#'
#' @param data_species 
#'
#' @return
#' @export
#'
#' @examples
plot_size_distribution_patch <- function(data_species) {
  
  data_species %>%
    dplyr::filter(!is.na(density)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = height, color = as.character(species),
                        group = cohort, alpha = density)) +
    ggplot2::geom_line() + 
    ggplot2::labs(x = "Time (years)",
         y = "Height (m)",
         color = "Species", 
         alpha = "Density (/m/m2)") +
    ggplot2::theme_classic()
}



#' Title
#'
#' @param patch_output_data 
#' @param steps 
#'
#' @return
#' @export
#'
#' @examples
plot_patch_output_across_time <- function(patch_output_data,
                                          steps = c(100, 200, 250, 300, 308)) 
  {
  patch_output_data %>% 
    dplyr::mutate(time_easy = format(time, digits = 0)) %>% 
    dplyr::filter(step %in% steps) %>% 
    ggplot2::ggplot(ggplot2::aes(x = canopy_openness, y = height, group = time_easy)) +
    ggplot2::geom_line() +
    ggplot2::theme_light() + 
    ggplot2::facet_wrap(~time_easy)
}

