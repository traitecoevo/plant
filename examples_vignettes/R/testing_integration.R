
expand_input_data <- function(...){    
  params = tibble(...)  
  
  params %>%
    mutate(b = calc_vul_b(p_50, c)) %>%
    mutate(psi_crit = calc_psi_crit(b,c)) %>%
    mutate(K_s = p_50_2_K_s(p_50)) %>%
    mutate(k_l_max = calc_k_l_max(K_s, huber_value, h)) %>%
    mutate(sun_rise = sun_rise(day, latitude),
           sun_down = sun_rise + (12-sun_rise)*2,
           day_length = sun_down - sun_rise,
           day_floor = floor(day),
           day_length_dec = day_length/24) 
}

create_mult_point <- function(...){
  data = tibble(...) 
  day_length = data$day_length
  integral_lengths = day_length/data$focal_points
  first_point = data$sun_rise + integral_lengths/2
  last_point = data$sun_down - integral_lengths/2
  
  if(data$focal_points == 1){
    return(first_point)
  } else{
    seq(first_point, last_point, by = integral_lengths)
  }
}

make_leaf <- function(...){
  params <- tibble(...)
  
  leaf <- plant:::Leaf(vcmax = params$vcmax, p_50 = params$p_50, c = params$c, b = params$b, psi_crit = params$psi_crit, K_s = params$K_s, epsilon_leaf = 0.0001, beta1 = params$beta1, beta2 = 1.5)
  
  params$sapwood_volume_per_area <- params$huber_value * params$h
  
  leaf$set_physiology(PPFD = params$PAR*params$E, psi_soil = params$psi_soil, leaf_specific_conductance_max = params$k_l_max, atm_vpd = params$VPD_hr, ca = params$ca, sapwood_volume_per_leaf_area = params$sapwood_volume_per_area)
  return(leaf)
}

leaf_states_rates_from_leaf<- function(leaf){
  leaf$optimise_psi_stem_Sperry()
  output_tibble <- tibble(profit = leaf$profit_, transpiration = leaf$transpiration_, psi_stem = leaf$opt_psi_stem_, g_c = leaf$stom_cond_CO2_)
}


visualise_integrals <- function(variable, ...){
  data_to_plot<- tibble(...)  
  

  if(variable == "profit"){

  data_to_plot  %>% unnest(coords) %>% unnest(coords) %>% ggplot() + geom_line(data = data_to_plot %>% filter(focal_points == 60) %>% select(-focal_points), aes(x=instant_through_day_hr , y = profit)) +
    geom_polygon(aes(x = x, y= profit_y), alpha = 0.3, col = "orange" , fill = "orange") +
    geom_point(aes(x = instant_through_day_hr, y= profit), col = "orange" , fill = "orange") +
    facet_wrap(~focal_points) +
    theme_classic() +
      geom_text(aes(x = 12, y= Inf, label = paste("n = ", focal_points)), vjust = 1, size =6)
  } else{
    data_to_plot  %>% unnest(coords) %>% unnest(coords) %>% ggplot() + geom_line(data = data_to_plot %>% filter(focal_points == 60) %>% select(-focal_points), aes(x=instant_through_day_hr , y = transpiration)) +
      geom_polygon(aes(x = x, y= transpiration_y), alpha = 0.3, col = "orange" , fill = "orange") +
      geom_point(aes(x = instant_through_day_hr, y= transpiration), col = "orange" , fill = "orange") +
      facet_wrap(~focal_points) +
      theme_classic()+
      geom_text(aes(x = 12, y= Inf, label = paste("n = ", focal_points)), vjust = 1, size =6)
  }
}

make_polygon_coordinates2 <- function(...){

  data <- tibble(...)
  
  data %>%
    select(transpiration, profit, focal_points, day_length, instant_through_day_hr, sun_rise, sun_down, day_length, run) %>%
    ungroup %>%
    mutate(coords = pmap(., create_polygon_coordinates)) %>%
    select(coords)
}

create_polygon_coordinates <- function(...){
  
  data2 <- tibble(...)
  day_length = data2$day_length
  integral_lengths = day_length/data2$focal_points
  
  middle_point = integral_lengths/2 + integral_lengths*(data2$run - 1) + data2$sun_rise
  lower_point = middle_point - integral_lengths/2
  higher_point = middle_point + integral_lengths/2
  
  data3 <- data2 %>% mutate(lower_point = lower_point,
                            higher_point = higher_point)
  
  data3 %>%
    tibble(x = c(.$lower_point, .$lower_point, .$higher_point, .$higher_point), profit_y = c(0, .$profit, .$profit, 0), transpiration_y = c(0, .$transpiration, .$transpiration, 0)) %>%
    select(x, profit_y, transpiration_y)
}





plotting_function <- function(..., variable){
  data <- tibble(...)
  col_name_continuous <- data$continuous_parameter[1]
  col_name_factorial <- data$factorial_parameter[1]
  
  col_name_factorial2 <- data$factorial_parameter2[1]
  col_name_continuous2 <- data$continuous_parameter2[1]
  
  if(variable == "profit"){
  
  ggplot(data)+
    geom_line(aes(x = .data[[col_name_continuous]], y = daily_profit_integrated, group = interaction(as.factor(.data[[col_name_factorial]]), focal_points), colour = as.factor(.data[[col_name_factorial]]), linetype = as.factor(focal_points)), show.legend = FALSE)+
    theme_classic() +
    ylab("")-> p
  }    
  
  if(variable == "transpiration"){
    
    ggplot(data)+
      geom_line(aes(x = .data[[col_name_continuous]], y = daily_transpiration_integrated, group = interaction(as.factor(.data[[col_name_factorial]]), focal_points), colour = as.factor(.data[[col_name_factorial]]), linetype = as.factor(focal_points)), show.legend = FALSE)+
      theme_classic() +
      ylab("") -> p

  }
  
  p + xlab(label = col_name_continuous2) -> p
  # p + scale_linetype_discrete(name = "Focal points") -> p
  
  if(!col_name_factorial %in% c("vcmax","p_50","huber_value")){
    # p + scale_colour_discrete(name = col_name_factorial2) -> p
  }
  
  if(col_name_continuous == "ca"){
    p + xlab(expression(C[a]~(Pa))) -> p
  }
  if(col_name_continuous == "psi_soil"){
    p + xlab(expression(psi[soil]))-> p
  }
  if(col_name_continuous == "temp_diff"){
    p + xlab(expression(Delta~T~(degree*C)))-> p
  }
  
  if(col_name_factorial == "vcmax"){
    # p + scale_colour_discrete(name = expression(V[c,max]~(mu*mol~m^{-2}~s^{-1}))) -> p
  }
  if(col_name_factorial == "p_50"){
    # p + scale_colour_discrete(name = expression(P[50]~(-MPa))) -> p
  }
  if(col_name_factorial == "huber_value"){
    # p + scale_colour_discrete(name = expression(HV~(m^2~Sapwood~area~m^{-2}~Leaf~area)))-> p
  }
  
  p+theme(text = element_text(size = 14))  
  p + theme(legend.position = "none")  + guides(colour = FALSE, linetype = FALSE)
  return(p)
}




