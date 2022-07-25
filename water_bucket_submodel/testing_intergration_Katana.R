
test2 <- 1
saveRDS(test2, "plant/water_bucket_submodel/test2.RDS")

rm(list=ls())

packages <- c("dplyr","tidyr","ggplot2","purrr")

install.packages(packages)

test3 <- 1
saveRDS(test3, "plant/water_bucket_submodel/test3.RDS")

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(deSolve)

test4 <- 1
saveRDS(test4, "plant/water_bucket_submodel/test4.RDS")


install.packages("devtools")
library(devtools)


test5 <- 1
saveRDS(test5, "plant/water_bucket_submodel/test5.RDS")


packages <- c("Rcpp", "R6", "crayon", "nleqslv", "BB" ,"BH")
install.packages(packages)


test6 <- 1
saveRDS(test6, "plant/water_bucket_submodel/test6.RDS")



devtools::install_deps("plant", upgrade = TRUE)

devtools::load_all("plant")

source("plant/water_bucket_submodel/R/solar_model.R")
source("plant/water_bucket_submodel/R/reference_leaf_updated.R")
source("plant/water_bucket_submodel/R/reference_soil.R")

test1 <- 1
saveRDS(test1, "plant/water_bucket_submodel/test1.RDS")

calc_VPH <- function(instant_of_day, VP_9, VP_15){
  gap = 18/24
  
  if(instant_of_day <= 9/24){
    VP_15 + (VP_9 - VP_15)*(9/24+instant_of_day)/gap  
  } else if(instant_of_day > 9/24 & instant_of_day <=15/24){
    VP_9 + (VP_15 - VP_9)*(instant_of_day - 9/24)/(15/24-9/24)
  } else if (instant_of_day > 15/24){
    VP_15 + (VP_9 - VP_15)*(instant_of_day - 15/24)/gap
  }
}

calc_temp <- function(instant_of_day, sun_up, sun_down, day_length, min_temp, max_temp){
  
  instant_of_day <- instant_of_day*24
  sun_up <- sun_up*24
  sun_down <- sun_down*24
  day_length <- day_length*24
  
  a_t = 1.86;
  b_t = 2.2;     #nighttime coeffcient
  c_t = -0.17;   #lag of the min temp from the time of runrise
  
  day_length = day_length 
  night_length = 24 - day_length
  
  m = sun_down - sun_up + c_t
  tset = (max_temp - min_temp) * sin(pi * m / (day_length + 2.0 * a_t)) + min_temp
  
  
  m = instant_of_day - sun_up + c_t
  
  if(instant_of_day >= sun_up & instant_of_day <= sun_down){
    min_temp + (max_temp - min_temp) * sin(pi * m / (day_length + 2.0 * a_t)) 
  }
  else{
    if(instant_of_day > sun_down){
      n = instant_of_day - sun_down
    }
    else if(instant_of_day < sun_up){
      n = 24 + instant_of_day - sun_down
    }
    d = (tset - min_temp)/(exp(b_t)-1)
    
    (min_temp - d) + (tset - min_temp - d) * exp(-b_t*n/(night_length + c_t))
  }
}
test8 <- 1
saveRDS(test8, "plant/water_bucket_submodel/test8.RDS")

calc_daily_model <- function(params){
  params %>%
    mutate(daily_data = map2(latitude, day, ~calc_PPFD_instant(.x, .y))) %>%
    select(-c(latitude, day)) %>%
    unnest(daily_data) %>%
    rowwise() %>%
    mutate(VPH = calc_VPH(instant_of_day, VPD_9, VPD_15), 
           temp = calc_temp(instant_of_day, sun_up, sun_down, day_length, min_temp = air_temp_c-diff_temp, max_temp = air_temp_c + diff_temp)) %>%
    mutate(VPD_hourly = max(0.001,calc_vpsat(temp) - VPH)) %>%
    mutate(profit = find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD_umol_m2_s_instant, b = b, c = c, psi_crit = psi_crit, atm_vpd = VPD_hourly, ca = ca)$profit_root)
}

calc_daily_model <- function(params){
  cluster <- new_cluster(2)
  params1 <- params %>% group_by(day) %>% partition(cluster)
  params1 %>%
    mutate(daily_data = map2(latitude, day, ~calc_PPFD_instant(.x, .y))) %>%
    select(-c(latitude, day)) %>%
    unnest(daily_data) %>%
    rowwise() %>%
    mutate(VPH = calc_VPH(instant_of_day, VPD_9, VPD_15), 
           temp = calc_temp(instant_of_day, sun_up, sun_down, day_length, min_temp = air_temp_c-diff_temp, max_temp = air_temp_c + diff_temp)) %>%
    mutate(VPD_hourly = max(0.001,calc_vpsat(temp) - VPH)) %>%
    mutate(profit = find_opt_psi_stem(psi_soil = psi_soil, k_l_max = k_l_max, vcmax = vcmax, PPFD = PPFD_umol_m2_s_instant, b = b, c = c, psi_crit = psi_crit, atm_vpd = VPD_hourly, ca = ca)$profit_root) %>%
    collect()
}


calc_one_point <- function(data){
  data %>%
    mutate(midday = (sun_down - sun_up)/2 + sun_up) %>%
    filter(abs(instant_of_day - midday) == min(abs(instant_of_day - midday))) %>%
    pull(profit) -> day_profit
  
  data %>%
    slice(1) %>%
    pull(profit)  -> night_profit
  day_night_profit <- c(day_profit, night_profit)
  
  data %>%
    slice(1) %>%
    summarise(day_length = day_length, night_length = 1 - day_length) -> day_night_length
  # sum(day_night_length*day_night_profit)
  return(sum(day_night_length*day_night_profit))
}

calc_three_point <- function(data){
  data %>%
    ungroup() %>%
    mutate(midday = (sun_down - sun_up)/2 + sun_up) %>%
    slice(which.min(abs(instant_of_day - midday))) -> day_profit
  
  data %>%
    ungroup() %>%
    mutate(after_sunrise = sun_up + 1/24) %>%
    slice(which.min(abs(instant_of_day - after_sunrise))) -> sunrise_profit
  
  data %>%
    ungroup() %>%
    mutate(before_sunset = sun_down - 1/24) %>%
    slice(which.min(abs(instant_of_day - before_sunset))) -> sunset_profit
  
  x <- c(sunrise_profit$sun_up, sunrise_profit$instant_of_day, day_profit$instant_of_day, sunset_profit$instant_of_day, sunset_profit$sun_down)
  
  y <- c(0, sunrise_profit$profit, day_profit$profit, sunset_profit$profit, 0)
  
  outputs %>%
    ungroup() %>%
    slice(1) %>%
    mutate(night_length = 1- day_length, night_profit = night_length*profit) %>%
    pull(night_profit) -> night_profit
  
  
  plant:::trapezium(x,y) + night_profit 
  
}

length_required <- 2

params <- set_params_soil(day = round(seq(1, 180, length.out = length_required)), latitude = seq(0, 60, length.out = length_required), vcmax = seq(30, 100, length.out = length_required), diff_temp = seq(5,10,length.out = length_required), VPD_9 = seq(0.5,2, length.out = length_required), diff_VH = seq(0.5,2, length.out = length_required), p50 = seq(1,5, length.out = length_required))

saveRDS(params, "plant/water_bucket_submodel/params.RDS")

outputs <- calc_daily_model(params)


saveRDS(outputs, "plant/water_bucket_submodel/outputs.RDS")

outputs %>%
  ungroup() %>%
  nest(daily_data = c(profit, instant_of_day, PPFD_umol_m2_s_instant, VPD_hourly, sun_up, sun_down, day_length, VPH, temp)) %>%
  mutate(true_daily_profit = map_dbl(daily_data, ~plant:::trapezium(.x$instant_of_day, .x$profit)),
         one_point_estimate = map_dbl(daily_data, calc_one_point),
         three_point_estimate = map_dbl(daily_data, calc_three_point)) -> outputs_integrated

saveRDS(outputs_integrated, "plant/water_bucket_submodel/outputs_integrated.RDS")
