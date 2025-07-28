

#' Title
#'
#' @param data_species ??
#'
#' @export
#' @importFrom rlang .data
plot_size_distribution <- function(data_species) {

  #relatavise a metric
  rel <- function(x, xmin = -4) {
    x[x < xmin] <- xmin
    xmax <- max(x, na.rm = TRUE)
    (x - xmin) / (xmax - xmin)
  }

  data_species %>%
    dplyr::filter(!is.na(.data$density)) %>%
    dplyr::mutate(relative_log_density = rel(.data$log_density)) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$height, color = as.character(.data$species),
                        group = interaction(.data$node, .data$species), alpha = .data$relative_log_density)) +
    ggplot2::geom_line() + 
    ggplot2::labs(x = "Time (years)",
         y = "Height (m)",
         color = "Species", 
         alpha = "Relative log(Density [/m/m2])") +
    ggplot2::scale_alpha(range=c(0.01,1))+ 
    ggplot2::theme_classic()
}


