E_from_Soil_to_Root_Collar <- function(P_x_r = P_x_r, P_soil,
                                       n_soil = 2, z_soil_mid = c(0.05, 0.15), dz = 0.1,
                                       LA = 1, c_r_H = c(20,30), c_r_V = c(20,30), beta_R_H = 3.4e3, beta_R_V = 9.4e4,
                                       rho, g){
  # Integration steps
  n = 20;
  
  c_r = c_r_H + c_r_V; #total root carbon in [mol C]
  r_R_H_min = beta_R_H / c_r_H #[MPa * s * (mol H2O)^-1]
  
  # Vertical root resistance in a given layer
  r_R_V = beta_R_V * ((dz^2)/ c_r_V)
  
  # Cumulative sum of vertical root resistance
  r_R_V_sum = cumsum(r_R_V)
  
  # Set up vector of root water uptake from layer
  E_soil = rep(NA, n_soil)
  
  # Set up vector of fractional resistance from each layer
  f_r = rep(NA, n_soil)
  
  for(i in 1:n_soil){
    
    # Find the most negative soil potential out of the given soil layer and the root collar
    P_src_min = min(P_soil[i], P_x_r);
    # Find the least negative soil potential out of the given soil layer and the root collar
    P_src_max = max(P_soil[i], P_x_r);
    
    # If root collar soil water potential equals the soil water potential in a given layer
    if(P_x_r == P_soil[i]){
      # Fraction of conductance in roots in a given layer at most negative soil water potential
      f_ri = VC_r(P_src_min)
      # Fraction of conductance in roots in a given layer at most negative soil water potential
      r_R_H = r_R_H_min[i] / f_ri; # [MPa * s * (mol H2O)^-1]
      # Total root resistance (horizantal plus vertical)
      r_R = r_R_H + r_R_V_sum[i]
      # Transpiration is equivalent to gravitational water loss (i.e. layer gains water)
      E_i = - (rho * g * z_soil_mid[i] / 10^6) / r_R / LA
    }
    else if((P_soil[i] - P_x_r) == (rho * g * z_soil_mid[i] / 1e6)){
      # If pressure difference perfectly balances gravity transpiration is equal to zero
      E_i = 0; # [mol H2O / m^2 / s]
    } 
    else{
      
      # Sequence through the most negative to least negative soil water potential
      P_src = seq(P_src_min, P_src_max, length.out = n);
      
      # Fraction of conductance in roots in a given layer at given soil water potential
      f_ri = VC_r(P_src)
      # Wherever the soil water potential is greater than zero, set the fractional loss of conductance to zero
      f_ri[P_src > 0] = VC_r(0);
      # Find the average f_ri
      f_ri = sum(f_ri) / n;
      
      # Find the horizantal resistance in a given layer by dividing the minimum resistance (i.e. maximum conductivity)
      # by the fractional loss of conductivity
      r_R_H = r_R_H_min[i] / f_ri; # [MPa * s * (mol H2O)^-1]
      
      # Find the total resistance in a given layer by adding the vertical resistance in that layer
      r_R = r_R_H + r_R_V_sum[i]; # [MPa * s * (mol H2O)^-1]
      
      # Transpiration is equal to the potentail gradient between the root collar and the soil, accounting for 
      # gravitational potential
      E_i = (P_soil[i] - P_x_r - rho * g * z_soil_mid[i]/ 10^6) / r_R / LA; # [mol H2O / m^2 / s]
    }
    
    # Add values to vector
    E_soil[i] = E_i
    f_r[i] = f_ri
  }
  
  # Total transpiration equal to sum of uptake from each layer
  E = sum(E_soil); # [mol H2O / m^2 / s]
  
  # Recalculate resistances in each layer (TODO: Bit unsure about why z_soil_mid is [i])
  r_R = (P_soil - P_x_r - rho * g * z_soil_mid / 10^6) / E_soil / LA; 
  r_R_H = r_R - r_R_V_sum; 
  r_R_H = pmax(r_R_H, r_R_H_min) #[MPa * s * (mol H2O)^-1]
  
  return(list(E, E_soil, r_R_H, r_R_V, f_r))
  }

# Vulnerability curve for leaf (parameterised from Potkay et al. 2021) [frac]
VC_l = function(psi){
  b_r = 0.85 #Mpa
  c_r = 0.81

  return(exp(-((-psi/b_r)^c_r)))
}
# Vulnerability curve for root (parameterised from Potkay et al. 2021) [frac]
VC_r = function(P_soil){
  b_r = 1.29 #Mpa
  c_r = 2.65
  return(exp(-((-P_soil/b_r)^c_r)))
}

# Vulnerability curve for sapwood (parameterised from Potkay et al. 2021) [frac]
VC_sw = function(psi){
  b_r = 5.32 #Mpa
  c_r = 0.80
  return(exp(-((-psi/b_r)^c_r)))
}
