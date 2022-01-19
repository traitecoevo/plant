

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
    mutate(relative_log_density = rel(log_density)) %>%
    ggplot2::ggplot(ggplot2::aes(x = time, y = height, color = as.character(species),
                        group = interaction(cohort,species), alpha = relative_log_density)) +
    ggplot2::geom_line() + 
    ggplot2::labs(x = "Time (years)",
         y = "Height (m)",
         color = "Species", 
         alpha = "Relative log(Density [/m/m2])") +
    ggplot2::scale_alpha(range=c(0.01,1))+ 
    ggplot2::theme_classic()
}

#' Relativise a score 
#'
#' @param x 
#' @param xmin 
#'
#' @return
#' @export
#'
#' @examples
rel <- function(x, xmin = -4) {
  x[x < xmin] <- xmin
  xmax <- max(x, na.rm=TRUE)
  (x - xmin) / (xmax - xmin)
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

