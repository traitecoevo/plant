x <- seq(0, 40, 1)

a <- -1
b <- 1
c <- 400

quadratic_equation <- function(x){
  a*(x-10)^2 + b*(x) + c
}

y <- quadratic_equation(x)

plot(y~x)
intial_ = 20
diff_initial = 0.1

compute_fst_deriv <- function(value, diff_value){
  x1 = value - diff_value
  x2 = value + diff_value
  
  f_x1 <- quadratic_equation(x1)
  f_x2 <- quadratic_equation(x2)
  
  rise = f_x2 - f_x1
  run = x2 - x1
  return(rise/run)
}

compute_second_deriv <- function(value, diff_value){
  x1 = value - diff_value
  x2 = value + diff_value
  x0 = value
  
  f_x0 <- quadratic_equation(x0)
  f_x1 <- quadratic_equation(x1)
  f_x2 <- quadratic_equation(x2)
  
  rise = f_x2 - 2*f_x0 + f_x1
  run = (x2 - x0)^2
  return(rise/run)
}

epsilon = 0.001
x_t1 = 20
x_t2 = x_t1 - compute_fst_deriv(x_t1, diff_initial)/compute_second_deriv(x_t1, diff_initial)
while(abs(x_t2 - x_t1) > epsilon){
  x_t2 = x_t1 - compute_fst_deriv(x_t1, diff_initial)/compute_second_deriv(x_t1, diff_initial)
  x_t1 = x_t2 - compute_fst_deriv(x_t2, diff_initial)/compute_second_deriv(x_t2, diff_initial)
  
}

plot(profits~seq(0, 3.54, length.out = 100))
points(x=2, l$calc_profit_Sperry(900, 0, 2, k_l_max), col = "red")
abline(a = 21.11455, b =-7.119663)

l$calc_profit_Sperry(900, 0, 2, k_l_max) - compute_fst_deriv(2, 0.01)*2

compute_fst_deriv <- function(value, diff_value, psi_soil){
  x1 = value - diff_value
  x2 = value + diff_value
  
  f_x1 <- l$calc_profit_Sperry(50, psi_soil, x1, k_l_max)
  f_x2 <- l$calc_profit_Sperry(50, psi_soil, x2, k_l_max)
  
  rise = f_x2 - f_x1
  run = 2*diff_value
  return(rise/run)
}

compute_second_deriv <- function(value, diff_value, psi_soil){
  x1 = value - diff_value
  x2 = value + diff_value
  x0 = value

  f_x0 <- l$calc_profit_Sperry(50, psi_soil, x0, k_l_max)
  f_x1 <- l$calc_profit_Sperry(50, psi_soil, x1, k_l_max)
  f_x2 <- l$calc_profit_Sperry(50, psi_soil, x2, k_l_max)

  rise = f_x2-2*f_x0+f_x1 
  run = diff_value^2
  return(rise/run)
}
diff_initial = 0.01
epsilon = 0.001

psi_soil = 2
x_t1 = 2.1
x_t1 = x_t2
x_t2 = x_t1 - compute_fst_deriv(x_t1, diff_initial, psi_soil)/compute_second_deriv(x_t1, diff_initial, psi_soil)
while(abs(x_t2 - x_t1) > epsilon){
  x_t1 = x_t2
  x_t2 = x_t1 - compute_fst_deriv(x_t1, diff_initial, psi_soil)/compute_second_deriv(x_t1, diff_initial, psi_soil)
  
}

plot
abline(a = 21.11455, b =-7.119663)

plot(map_dbl(seq(0, 3.54, length.out = 100), ~l$calc_profit_Sperry(50, 0, .x, k_l_max))~seq(0, 3.54, length.out = 100))
plot(map_dbl(seq(0, 3.54, length.out = 100), ~compute_fst_deriv(.x, 0.01))~seq(0, 3.54, length.out = 100))
points(x=2, compute_fst_deriv(2, 0.01), col = "red")

compute_second_deriv(2, 0.01)
compute_second_deriv(2, 0.01)
