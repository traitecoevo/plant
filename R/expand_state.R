expand_state <- function(results) {
  
  data <- 
    results$species %>% 
    split(., .$species)
  
  for(i in seq_len(results$n_spp)) {

    s <- results$p$strategies[[i]]
    s$eta_c <- 1 - 2/(1 + s$eta) + 1/(1 + 2*s$eta)
    
    data[[i]] <- 
      data[[i]] %>%
      mutate(
        # These are formulas from ff16_strategy.cpp 
        # ideally wouldn't have to copy them here
        # could we expose them from startegy object
        # and call them directly?
        
        area_leaf = (height / s$a_l1)^(1.0 / s$a_l2),
        mass_leaf = area_leaf * s$lma,
        area_sapwood = area_leaf * s$theta,
        mass_sapwood = area_sapwood * height * s$eta_c * s$rho,
        area_bark = s$a_b1 * area_leaf * s$theta,
        mass_bark = area_bark * height * s$eta_c * s$rho,
        area_stem = area_bark + area_sapwood + area_heartwood,
        diameter_stem = sqrt(4 * area_stem / pi),
        mass_root = s$a_r1 * area_leaf,
        mass_live = mass_leaf + mass_sapwood + mass_bark + mass_root,
        mass_total =  mass_leaf + mass_bark + mass_sapwood +  mass_heartwood + mass_root,
        mass_above_ground = mass_leaf + mass_bark + mass_sapwood +  mass_heartwood
      )
    }
  
  results$species <- data %>% bind_rows()
  
  results
}