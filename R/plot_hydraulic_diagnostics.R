# Diagnostic plot function for TF24 SCM results
# Creates all panels, assembles them with patchwork, and saves to file.

plot_diagnostics <- function(results, x, y,
                             filename = "diagnostic_plot.png",
                             height = 20, width = 15) {
  library(tidyverse)
  library(patchwork)

  # Derive total patch leaf area over time (used in transpiration_per_leaf_area)
  total_area_leaf <- results %>%
    expand_state() %>%
    pluck("species") %>%
    integrate_over_size_distribution() %>%
    pull(area_leaf)

  # Patch leaf area over time
  results %>%
    expand_state() %>%
    pluck("species") %>%
    integrate_over_size_distribution() %>%
    ggplot(aes(time, area_leaf)) +
    geom_line(linewidth = 2, colour = "forestgreen", na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(Patch ~ leaf ~ area, " (", m^2, " ", m^{-2}, ")"))) +
    xlab("Time (yr)") -> patch_leaf_area

  # Size distribution over time
  results$species %>%
    drop_na() %>%
    plot_size_distribution() +
    theme_bw() +
    theme(text = element_text(size = 20)) -> size_distribution

  # Soil water potential by depth over time
  results$env$soil_moist %>%
    mutate(soil_depth = results$env$soil_depth$soil_depth) %>%
    mutate(psi_soil = (1.78e3 * (soil_moist / 0.428)^-6.57) / 1e6) %>%
    ggplot(aes(x = time, y = psi_soil, group = soil_depth, colour = soil_depth)) +
    geom_point(na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(psi[soil], " (", -Mpa, ")"))) +
    xlab("Time (years)") +
    labs(colour = "Soil depth (m)") -> soil_water_potential

  # Soil moisture by depth over time
  results$env$soil_moist %>%
    mutate(soil_depth = results$env$soil_depth$soil_depth) %>%
    mutate(psi_soil = (1.78e3 * (soil_moist / 0.428)^-6.57) / 1e6) %>%
    ggplot(aes(x = time, y = soil_moist, group = soil_depth, colour = soil_depth)) +
    geom_point(na.rm = TRUE) +
    geom_line(na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(psi[soil], " (", -Mpa, ")"))) +
    xlab("Time (years)") +
    labs(colour = "Soil depth (m)") -> soil_moisture

  # Stem water potential by height over time
  results$species %>%
    ggplot(aes(x = time, y = opt_psi_stem)) +
    geom_line(aes(colour = height, group = node), na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(psi[stem], " (", -Mpa, ")"))) +
    xlab("Time (years)") +
    labs(colour = "Height (m)") -> stem_water_potential

  # Root water potential by height over time
  results$species %>%
    ggplot(aes(x = time, y = opt_root_psi)) +
    geom_line(aes(colour = height, group = node), na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(psi[root], " (", -Mpa, ")"))) +
    xlab("Time (years)") +
    labs(colour = "Height (m)") -> root_water_potential

  # Profit by height over time
  results$species %>%
    ggplot(aes(x = time, y = profit)) +
    geom_line(aes(colour = height, group = node), na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(Profit, " (", mu, mol ~ m^{-2} ~ s^{-1}, ")"))) +
    xlab("Time (years)") +
    labs(colour = "Height (m)") -> profit

  # Stomatal conductance by height over time
  results$species %>%
    ggplot(aes(x = time, y = stom_cond_CO2)) +
    geom_line(aes(colour = height, group = node, alpha = density), na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(g[s], " (", mol ~ m^{-2} ~ s^{-1}, ")"))) +
    xlab("Time (years)") +
    labs(colour = "Height (m)") -> stomatal_conductance

  # Rainfall time series
  tibble(x = x, y = y) %>%
    ggplot(aes(x = x, y = y)) +
    geom_line(na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(Rainfall, " (", m ~ yr^{-1}, ")"))) +
    xlab("Time (years)") -> rainfall

  # Transpiration over time
  tibble(
    deplet = results$env$soil_moist_cumulative_flux$sum_resource_depletion,
    time   = results$env$soil_moist_cumulative_flux$time
  ) %>%
    # Increment of the *cumulative* depletion divided by the interval width:
    # without /dt this is "depletion per output interval", which scales with the
    # (non-uniform) node spacing and looks janky. Dividing by dt gives the true
    # rate the y-axis claims (m yr^-1). See traitecoevo/plant#474.
    mutate(deplet = (deplet - lag(deplet)) / (time - lag(time))) %>%
    ggplot(aes(x = time, y = deplet)) +
    geom_line(na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(Transpiration, " (", m ~ yr^{-1}, ")"))) +
    xlab("Time (years)") -> transpiration

  # Transpiration per leaf area over time
  tibble(
    deplet = results$env$soil_moist_cumulative_flux$sum_resource_depletion,
    time   = results$env$soil_moist_cumulative_flux$time
  ) %>%
    # Increment of the *cumulative* depletion divided by the interval width:
    # without /dt this is "depletion per output interval", which scales with the
    # (non-uniform) node spacing and looks janky. Dividing by dt gives the true
    # rate the y-axis claims (m yr^-1). See traitecoevo/plant#474.
    mutate(deplet = (deplet - lag(deplet)) / (time - lag(time))) %>%
    mutate(deplet_per_area = deplet / c(NA, total_area_leaf)) %>%
    ggplot(aes(x = time, y = deplet_per_area)) +
    geom_line(na.rm = TRUE) +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    ylab(expression(paste(Transpiration ~ per ~ leaf ~ area, " (", m ~ m^{-2} ~ yr^{-1}, ")"))) +
    xlab("Time (years)") -> transpiration_per_leaf_area

  # Assemble with patchwork and save
  combined <- (patch_leaf_area + size_distribution) /
    ((soil_water_potential + root_water_potential) / (stem_water_potential + profit)) /
    (stomatal_conductance + transpiration) /
    (transpiration_per_leaf_area + rainfall)

  ggsave(filename, plot = combined, height = height, width = width)

  print(combined)
}
