library(GPfit)

# Source: https://github.com/bearloga/bayesopt-tutorial-r/blob/master/BayesOpt.R

bayesian_optimize <- function(fn, 
                              evals = NULL,
                              iter = 1,
                              n_seed_evals = 4,
                              n_eval_pts = 1000,
                              minimize = TRUE, 
                              control = NULL, 
                              bounds = c(0, 1)) {
  
  
  if (is.null(control)) {
    control <- list(cov = list(type = "exponential", power = 1.95))
    # control <- list(cov = list(type = "matern", nu = 5/2))
  }
  
  # need seed evaluations to fit GP
  if(is.null(evals)) {
    evals <- data.frame(x = seq(bounds[1], bounds[2], 
                                length.out = n_seed_evals))
    
    evals$y <- map_dbl(evals$x, calibrate)
  }
 
  # TODO: GPfit wants x-scale on (0, 1]
  x_new = seq(bounds[1], bounds[2], len = n_eval_pts)

  for(i in 0:iter) {
    
    # current best eval
    if (minimize) {
      y_best <- min(evals$y, na.rm = TRUE)
    } else {
      y_best <- max(evals$y, na.rm = TRUE)
    }
    
    # evaluate next place to sample
    if(i > 0) {
      # Expected improvement acquisition function
      acq <- purrr::map2_dbl(mu, sigma, function(m, s) {
        if (s == 0) return(0)
        gamma <- (y_best - m) / s
        if (minimize) {
          phi <- pnorm(gamma)
        } else {
          phi <- 1 - pnorm(gamma)
        }
        return(s * (gamma * phi + dnorm(gamma)))
      })
      
      x_next <- x_new[which.max(acq)]
      y_next <- fn(x_next)
      
      evals <- rbind(evals, c(x_next, y_next))
    }
    
    # fit/re-fit emulator
    fit <- GP_fit(
      X = evals$x,
      Y = evals$y,
      corr = control$cov
    )
    
    preds <- predict.GP(fit, xnew = data.frame(x = x_new))
    mu <- preds$Y_hat
    sigma <- sqrt(preds$MSE)
  }
  
  minima = x_new[which.min(preds$Y_hat)]
  
  plot(preds$complete_data, type = "l", 
       main = paste("Expected improvement evaluations;",
                    "n = ", nrow(evals)),
       xlab = "B_lf1", ylab = "Error")
  
  points(evals, pch = 16)
  
  abline(v = minima, col = "red")
  legend("topleft", legend = "optima",
         pch = "|", col = "red")
  
  return(list(fit = fit, 
              evaluations = evals,
              minima = minima))
}





