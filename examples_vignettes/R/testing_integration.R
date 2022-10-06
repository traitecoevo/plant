
make_leaf <- function(...){
  params <- tibble(...)
  leaf <- plant:::Leaf(vcmax = params$vcmax, p_50 = params$p_50, c = params$c, b = params$b, psi_crit = params$psi_crit, huber_value = params$huber_value, K_s = params$K_s, epsilon_leaf = 0.0001)
  leaf$set_physiology(PPFD = params$PAR*params$E, psi_soil = params$psi_soil, k_l_max = params$k_l_max, atm_vpd = params$VPD_hr, ca = params$ca)
  return(leaf)
}


leaf_states_rates_from_leaf<- function(leaf){
  leaf$optimise_psi_stem_Sperry()
  output_tibble <- tibble(profit = leaf$profit, GSS_count = leaf$GSS_count, count = leaf$count, method = leaf$method, transpiration = leaf$E, psi_stem = leaf$opt_psi_stem, g_c = leaf$g_c)
}

calc_one_point <- function(data){
  (12 - min(data$instant_through_day_hr)) * 2 -> day_length_hr
  data %>%
    filter(abs(instant_through_day_dec - 0.5) == min(abs(instant_through_day_dec - 0.5))) %>%
    mutate(daily_profit = profit*(day_length_hr)) %>%
    sample_n(1) -> out
  
  out$daily_profit
}

extract_time_slice <- function(quantile, data){
  
  data %>%
    filter(abs(outputs$instant_through_day_hr - quantile(outputs$instant_through_day_hr, quantile)) %in% min(abs(outputs$instant_through_day_hr - quantile(outputs$instant_through_day_hr, quantile)))) %>%
    sample_n(1)
}


calc_mult_point <- function(...){
  data = tibble(...) %>%
    nest(outputs = outputs)
  first_num = 1/data$focal_points/2 
  quantiles = seq(from = first_num, 1, by = first_num*2)
  
  data %>%
    unnest(outputs) -> data
  
  map(quantiles, extract_time_slice, data) %>%
    bind_rows() %>%
    mutate(profit_integrated = outputs$profit * day_length/focal_points,
           transpiration_integrated = outputs$transpiration * day_length/focal_points) %>% 
    unnest(outputs) %>% 
    select(profit, profit_integrated, transpiration, transpiration_integrated, instant_through_day_hr) %>% 
    rename(transpiration_point = transpiration, profit_point = profit, instant_through_day_hr_point = instant_through_day_hr) %>% 
    nest(cols = c(transpiration_point, profit_point, instant_through_day_hr_point))
}



profit_curves_from_leaf<- function(leaf, psi_crit){
  tibble(inst_psi_stem = seq(leaf$psi_soil_, psi_crit, length.out = 100)) %>%
    rowwise() %>%
    mutate(inst_profit = leaf$profit_psi_stem_Sperry(inst_psi_stem),
           inst_cost = leaf$lambda_ * leaf$calc_hydraulic_cost_Sperry(inst_psi_stem),
           inst_ben = inst_profit + inst_cost)
}


visualise_integrals <- function(variable, ...){
  data_to_plot<- tibble(...)  
  
  if(variable == "profit"){
  data_to_plot  %>% unnest(coords) %>% unnest(coords) %>% unnest(outputs) %>% ggplot(aes(x=instant_through_day_hr , y = profit)) + geom_line() +
    geom_polygon(aes(x = x, y= profit_y), alpha = 0.3, col = "orange" , fill = "orange") +
    geom_point(aes(x = instant_through_day_hr_point, y= profit_point), col = "orange" , fill = "orange") +
    facet_wrap(~focal_points) +
    theme_classic()
  } else{
  data_to_plot  %>% unnest(coords) %>% unnest(coords) %>% unnest(outputs) %>% filter(focal_points == 50) %>% ggplot(aes(x=instant_through_day_hr , y = transpiration)) + geom_line() +
    geom_polygon(aes(x = x, y= transpiration_y), alpha = 0.3, col = "orange" , fill = "orange") +
    geom_point(aes(x = instant_through_day_hr_point, y= transpiration_point), col = "orange" , fill = "orange") +
    facet_wrap(~focal_points) +
    theme_classic()
  }
}




make_polygon_coordinates <- function(...){
  data <- tibble(...)
  
  
  data %>%
    nest(outputs = outputs) %>%
    select(transpiration_point, profit_point, focal_points, day_length, instant_through_day_hr_point, sun_rise, sun_down, day_length, run) %>%
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
    tibble(x = c(.$lower_point, .$lower_point, .$higher_point, .$higher_point), profit_y = c(0, .$profit_point, .$profit_point, 0), transpiration_y = c(0, .$transpiration_point, .$transpiration_point, 0)) %>%
    select(x, profit_y, transpiration_y)
}

expand_input_data <- function(...){    
  params = tibble(...)  
  
  params %>%
    # mutate(c = 2.04) %>%
    mutate(b = calc_vul_b(p_50, c)) %>%
    mutate(psi_crit = calc_psi_crit(b,c)) %>%
    mutate(K_s = p_50_2_K_s(p_50)) %>%
    mutate(k_l_max = calc_k_l_max(K_s, huber_value, h)) %>%
    mutate(sun_rise = sun_rise(day, latitude),
           sun_down = sun_rise + (12-sun_rise)*2,
           day_length = sun_down - sun_rise,
           day_floor = floor(day),
           day_length_dec = day_length/24) 
  # %>%
    # expand_grid(instant_through_day_hr = seq(sun_rise, sun_down, length.out = 300)) %>%
    # mutate(instant_through_day_dec = instant_through_day_hr/24,
    #        dec_day_time = day_floor+ instant_through_day_dec) %>% 
    # mutate(solar_angle = solar_angle(dec_day_time, latitude),
    #        PAR = PAR_given_solar_angle(solar_angle)) %>%
    # rowwise() %>%
    # mutate(temp = calc_temp(instant_through_day_hr, day_length, mean_temp, temp_diff, sun_rise, sun_down)) %>%
    # mutate(VPD_hr = max(0.05, calc_vpsat(temp) - VP))
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

calc_mult_point2 <- function(...){
  data = tibble(...) %>%
    nest(outputs = outputs)
  first_num = 1/data$focal_points/2 
  quantiles = seq(from = first_num, 1, by = first_num*2)
  
  data %>%
    unnest(outputs) -> data
  
  map(quantiles, extract_time_slice, data) %>%
    bind_rows() %>%
    mutate(profit_integrated = outputs$profit * day_length/focal_points,
           transpiration_integrated = outputs$transpiration * day_length/focal_points) %>% 
    unnest(outputs) %>% 
    select(profit, profit_integrated, transpiration, transpiration_integrated, instant_through_day_hr) %>% 
    rename(transpiration_point = transpiration, profit_point = profit, instant_through_day_hr_point = instant_through_day_hr) %>% 
    nest(cols = c(transpiration_point, profit_point, instant_through_day_hr_point))
}

