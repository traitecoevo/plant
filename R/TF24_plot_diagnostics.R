#' Diagnostic panel plot for a TF24 stand
#'
#' Assembles a multi-panel diagnostic figure from a collected TF24 SCM run:
#' patch leaf area, the size distribution, soil water potential, stem and root
#' water potentials, profit, stomatal conductance, transpiration (total and per
#' leaf area) and the rainfall driver. Panels are combined with `patchwork`.
#'
#' `ggplot2` and `patchwork` are Suggested packages; the function errors if
#' either is unavailable. The returned plot is neither drawn nor saved — print
#' it or pass it to [ggplot2::ggsave()] as needed.
#'
#' @param results A collected TF24 SCM run, i.e. the output of
#'   `run_scm(..., collect = TRUE)` for a TF24 strategy/environment. Carries the
#'   tidied patch state (`results$species`) and the environment series
#'   (`results$env$soil_moist`, `soil_depth`, `soil_moist_cumulative_flux`).
#' @param x,y Numeric vectors giving the rainfall time series — time (years) and
#'   rainfall rate (m yr^-1) — as supplied to the environment's extrinsic
#'   "rainfall" driver.
#'
#' @return The assembled `patchwork` plot.
#'
#' @export
#' @importFrom rlang .data
TF24_plot_diagnostics <- function(results, x, y) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("plot_diagnostics() requires the 'ggplot2' and 'patchwork' packages.")
  }

  # Derive total patch leaf area over time (used in transpiration_per_leaf_area)
  total_area_leaf <- results %>%
    expand_state() %>%
    purrr::pluck("species") %>%
    integrate_over_size_distribution() %>%
    dplyr::pull(.data$area_leaf)

  # Patch leaf area over time
  results %>%
    expand_state() %>%
    purrr::pluck("species") %>%
    integrate_over_size_distribution() %>%
    ggplot2::ggplot(ggplot2::aes(.data$time, .data$area_leaf)) +
    ggplot2::geom_line(linewidth = 2, colour = "forestgreen", na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(Patch ~ leaf ~ area, " (", m^2, " ", m^{-2}, ")"))) +
    ggplot2::xlab("Time (yr)") -> patch_leaf_area

  # Size distribution over time
  results$species %>%
    tidyr::drop_na() %>%
    plot_size_distribution() +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) -> size_distribution

  # Soil water potential by depth over time
  results$env$soil_moist %>%
    dplyr::mutate(soil_depth = results$env$soil_depth$soil_depth) %>%
    dplyr::mutate(psi_soil = (1.78e3 * (.data$soil_moist / 0.428)^-6.57) / 1e6) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$psi_soil, group = .data$soil_depth, colour = .data$soil_depth)) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(psi[soil], " (", -Mpa, ")"))) +
    ggplot2::xlab("Time (years)") +
    ggplot2::labs(colour = "Soil depth (m)") -> soil_water_potential

  # Soil moisture by depth over time
  results$env$soil_moist %>%
    dplyr::mutate(soil_depth = results$env$soil_depth$soil_depth) %>%
    dplyr::mutate(psi_soil = (1.78e3 * (.data$soil_moist / 0.428)^-6.57) / 1e6) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$soil_moist, group = .data$soil_depth, colour = .data$soil_depth)) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(psi[soil], " (", -Mpa, ")"))) +
    ggplot2::xlab("Time (years)") +
    ggplot2::labs(colour = "Soil depth (m)") -> soil_moisture

  # Stem water potential by height over time
  results$species %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$opt_psi_stem)) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$height, group = .data$node), na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(psi[stem], " (", -Mpa, ")"))) +
    ggplot2::xlab("Time (years)") +
    ggplot2::labs(colour = "Height (m)") -> stem_water_potential

  # Root water potential by height over time
  results$species %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$opt_root_psi)) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$height, group = .data$node), na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(psi[root], " (", -Mpa, ")"))) +
    ggplot2::xlab("Time (years)") +
    ggplot2::labs(colour = "Height (m)") -> root_water_potential

  # Profit by height over time
  results$species %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$profit)) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$height, group = .data$node), na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(Profit, " (", mu, mol ~ m^{-2} ~ s^{-1}, ")"))) +
    ggplot2::xlab("Time (years)") +
    ggplot2::labs(colour = "Height (m)") -> profit

  # Stomatal conductance by height over time
  results$species %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$stom_cond_CO2)) +
    ggplot2::geom_line(ggplot2::aes(colour = .data$height, group = .data$node, alpha = .data$density), na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(g[s], " (", mol ~ m^{-2} ~ s^{-1}, ")"))) +
    ggplot2::xlab("Time (years)") +
    ggplot2::labs(colour = "Height (m)") -> stomatal_conductance

  # Rainfall time series
  tibble::tibble(x = x, y = y) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(Rainfall, " (", m ~ yr^{-1}, ")"))) +
    ggplot2::xlab("Time (years)") -> rainfall

  # Transpiration over time
  tibble::tibble(
    deplet = results$env$soil_moist_cumulative_flux$sum_resource_depletion,
    time   = results$env$soil_moist_cumulative_flux$time
  ) %>%
    # Increment of the *cumulative* depletion divided by the interval width:
    # without /dt this is "depletion per output interval", which scales with the
    # (non-uniform) node spacing and looks janky. Dividing by dt gives the true
    # rate the y-axis claims (m yr^-1). See traitecoevo/plant#474.
    dplyr::mutate(deplet = (.data$deplet - dplyr::lag(.data$deplet)) / (.data$time - dplyr::lag(.data$time))) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$deplet)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(Transpiration, " (", m ~ yr^{-1}, ")"))) +
    ggplot2::xlab("Time (years)") -> transpiration

  # Transpiration per leaf area over time
  tibble::tibble(
    deplet = results$env$soil_moist_cumulative_flux$sum_resource_depletion,
    time   = results$env$soil_moist_cumulative_flux$time
  ) %>%
    # Increment of the *cumulative* depletion divided by the interval width:
    # without /dt this is "depletion per output interval", which scales with the
    # (non-uniform) node spacing and looks janky. Dividing by dt gives the true
    # rate the y-axis claims (m yr^-1). See traitecoevo/plant#474.
    dplyr::mutate(deplet = (.data$deplet - dplyr::lag(.data$deplet)) / (.data$time - dplyr::lag(.data$time))) %>%
    dplyr::mutate(deplet_per_area = .data$deplet / c(NA, total_area_leaf)) %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$time, y = .data$deplet_per_area)) +
    ggplot2::geom_line(na.rm = TRUE) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 20)) +
    ggplot2::ylab(expression(paste(Transpiration ~ per ~ leaf ~ area, " (", m ~ m^{-2} ~ yr^{-1}, ")"))) +
    ggplot2::xlab("Time (years)") -> transpiration_per_leaf_area

  # Assemble with patchwork and save
  combined <- (patch_leaf_area + size_distribution) /
    ((soil_water_potential + root_water_potential) / (stem_water_potential + profit)) /
    (stomatal_conductance + transpiration) /
    (transpiration_per_leaf_area + rainfall)

  combined
}
