
#' Title
#'
#' @param results ??
#'
#' @return ???
#' @export
#' @importFrom rlang .data
expand_state <- function(results) {
  
  data <- split(results$species, results$species$species)
  
  for(i in seq_len(results$n_spp)) {

    s <- results$p$strategies[[i]]
    s$eta_c <- 1 - 2/(1 + s$eta) + 1/(1 + 2*s$eta)
    
    data[[i]] <- 
      data[[i]] %>%
      dplyr::mutate(
        # These are formulas from ff16_strategy.cpp 
        # ideally wouldn't have to copy them here
        # could we expose them from startegy object
        # and call them directly?
        
        area_leaf = (.data$height / s$a_l1)^(1.0 / s$a_l2),
        mass_leaf = .data$area_leaf * s$lma,
        area_sapwood = .data$area_leaf * s$theta,
        mass_sapwood = .data$area_sapwood * .data$height * s$eta_c * s$rho,
        area_bark = s$a_b1 * .data$area_leaf * s$theta,
        mass_bark = .data$area_bark * .data$height * s$eta_c * s$rho,
        area_stem = .data$area_bark + .data$area_sapwood + .data$area_heartwood,
        diameter_stem = sqrt(4 * .data$area_stem / pi),
        mass_root = s$a_r1 * .data$area_leaf,
        mass_live = .data$mass_leaf + .data$mass_sapwood + .data$mass_bark + .data$mass_root,
        mass_total =  .data$mass_leaf + .data$mass_bark + .data$mass_sapwood +  .data$mass_heartwood + .data$mass_root,
        mass_above_ground = .data$mass_leaf + .data$mass_bark + .data$mass_sapwood +  .data$mass_heartwood
      )
    }
  
  results$species <- data %>% dplyr::bind_rows()
  
  results
}
