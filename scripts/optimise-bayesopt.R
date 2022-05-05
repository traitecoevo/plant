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
      x = evals$x,
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


# source: https://bookdown.org/rbg/surrogates/chap7.html
library(laGP)
library(lhs)

eps <- sqrt(.Machine$double.eps) ## used lots below

# expected improvement of sampling x compared to the currently observed minima
EI <- function(gp, x, fmin) {
  if(is.null(nrow(x))) 
    x <- matrix(x, nrow=1)
  
  p <- predGPsep(gp, x, lite=TRUE)
  
  d <- fmin - p$mean
  sigma <- sqrt(p$s2)
  dn <- d / sigma
  ei <- d * pnorm(dn) + sigma * dnorm(dn)

  return(ei)
}


EI.search <- function(x, y, gp, bounds, multi_start = 5, tol = eps) {
  
  # find minima of data and start there
  m <- which.min(y)
  fmin <- y[m]
  start <- matrix(x[m, ], nrow=1)
  
  # randomly sample around starting point
  if(multi_start > 1) {
    start <- rbind(start, randomLHS(multi_start - 1, ncol(x)))
  } else {
    message("Using a multi-start optimisation routine
            can improve numerical instability")
  }
  
  max_ei = 0
  
  # repeatedly optimise
  for(i in 1:nrow(start)) {
    
    if(EI(gp, start[i, ], fmin) <= tol) {
      out <- list(value = -Inf)
      next 
    }
    
    # optim likes to minimise things, in this case the negative expected 
    # improvement of the next acquisition to improve the predicted minima of GP
    objective <- function(x, fmin, gp) {
      - EI(gp, x, fmin)
    }
    
    limits = tibble(bounds) %>% unnest_wider(bounds)
    
    out <- optim(start[i, ], objective, method = "L-BFGS-B", 
                 lower = limits$min,
                 upper = limits$max,
                 gp = gp, fmin = fmin)
    
    ei = -out$value

    # save best acquisition location
    if(ei > max_ei) {
      max_ei = ei
      x_next = matrix(out$par, nrow = 1)
    } 
  }

  return(x_next)
}

# @param bounds a list with upper and lower limits for each variable
generate_seeds <- function(f, bounds, n = 1,
                           target = NULL,
                           parameters = NULL,
                           evaluate = T) {
  x <- randomLHS(n, length(bounds))
  
  for(i in 1:length(bounds)) {
    x[, i] = bounds[[i]]$min + x[, i] * (bounds[[i]]$max - bounds[[i]]$min)
  }
  
  colnames(x) <- names(bounds)
  
  if(evaluate == T) {
    if(is.null(target) | is.null(parameters)) {
      stop("Target and parameters must be provided to evaluate seeds")
    }
    
    y <- apply(x, 1, f, target = target, parameters = parameters)
    return(list(x = x, y = y))
    
  } else {
    return(x)
  }
}


generate_range <- function(bounds, n = 1000) {
  x = matrix(nrow = n, ncol= length(bounds))
  basis = seq(0, 1, len = n)
  
  for(i in 1:length(bounds)) {
    x[, i] = bounds[[i]]$min + basis * (bounds[[i]]$max - bounds[[i]]$min)
  }
  
  return(x)
}


# @param x optimisable parameters
# @param target a function or vector of points to calculate the marginal likelihood 
# @param latitude used to specify site characteristics
# @param n_control_points number of points to evaluate LL if target is a function
calibrate <- function(x, target, latitude = 28.182,
                      parameters, 
                      environment = NULL,
                      n_control_pts = 10,
                      match_tissue = "mass_heartwood") {
  
  if(is.null(parameters$cohort_schedule_times)) {
    stop("Please run `create_schedule` first")
  }
  
  # assume optimising B_lf1 and hmat
  with(as.list(x), {
    site <- create_site(B_lf1 = B_lf1, latitude = latitude)
    
    sp <- create_species(hmat = hmat, site = site)
    # sp <- create_species(site = site)
    
    # overwrite species
    parameters$strategies <- sp$strategies
    parameters$birth_rate <- sp$birth_rate
    
    if(is.null(environment)) {
      environment = create_environment()
    }
    
    biomass = calculate_net_biomass(parameters, environment) 
    
    if(is.function(target)) {

      # drop final timepoint
      s = parameters$cohort_schedule_times[[1]]
      s = s[-length(s)]

      if(n_control_pts == 1) {
         # target endpoint
         pts <- max(s)
      } else { 
        # generate equidistant control points - not actually n :/
        pts <- s[findInterval(seq(0, max(s), len = n_control_pts), s)] %>%
          unique()
      }
      
      if(match_tissue == "total") {
        biomass <- group_by(biomass, time, species) %>%
          summarise(tissue = "total",
                    value = sum(value), 
                    .groups = "drop")
      }
      
      totals <- filter(biomass, 
                       tissue == match_tissue, 
                       time %in% pts) %>%
        mutate(target = target(time)) %>%
        mutate(error = target - value)
    } else {
      stop("Calibration for dataframes not yet implemented")
    }
    
    # loss fn - emphasise fit of total biomass by downweighting early ctrl pts
    weights = seq(0.1, 1, len = n_distinct(totals$time))
    error <- matrix(log(sum(abs(totals$error) * weights)), ncol = 1)
    
    return(error)
  })
}



optim.EI <- function(f, target, bounds, parameters, n_iter = 1, n_init = 6,
                     gp = NULL, evals = NULL) {

  if(is.null(bounds))
    stop("Please specify the bounds of optimisable parameters")
    
  if(is.null(evals)) {
    if(n_init < 1)
      stop("Provide a set of initial simulations (seed) or 
           randomly generate n_init new simulations")
    
    evals <- generate_seeds(f, n_init, target, bounds, parameters)
  }
  
  x <- evals$x
  y <- evals$y
  
  # initialise GP in C, with separable correlation for >1 dimensions
  # d = length scale, g = nugget, dK = save derivative for optimisation of d
  if(is.null(gp)) {
    gp_ptr <- newGPsep(x, y, d=0.1, g=1e-6, dK=TRUE)
    
    # priors for GP lengthscale, optimised using MLE
    da <- darg(list(mle = T), 
               generate_seeds(calibrate, bounds, n = 1000, evaluate = F))

  } else if(is.list(gp)) {
    gp_ptr = gp$pointer
    da <- gp$priors
  } else {
    stop("GP should have a pointer and a priors object")
  }
  
  ga <- garg(list(mle = T), y)
  
  mleGPsep(gp_ptr, param="d", tmin= da$min, tmax= da$max)  

  ## optimization loop of sequential acquisitions
  for(i in 1:n_iter) {
    
    # find best expected improvement
    x_next <- EI.search(x, y, gp_ptr, bounds)
    colnames(x_next) <- names(bounds)

    # collect data at proposed location
    y_next <- apply(x_next, 1, f, target = target, parameters = parameters)
    
    # update GP and length-scale
    updateGPsep(gp_ptr, x_next, y_next)
    mleGPsep(gp_ptr, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)
    
    # save new simulation
    x <- rbind(x, x_next)
    y <- c(y, y_next)
  }
  

  return(list(gp = list(pointer = gp_ptr, priors = da), 
              evaluations = list(x = x, y = y),
              minima = min(y)))
}


plot_preds <- function(gp_ptr, bounds, evals, par = 1, `3d` = FALSE) {
  # initial fit
  xx <- generate_range(bounds, n = 100)
  
  if(`3d`) {
    XX <- expand.grid(xx[, 1], xx[, 2])
    
    pred <- predGPsep(gp_ptr, XX, lite=TRUE)
    
    z <- matrix(pred$mean, ncol = nrow(xx))
    
    persp(xx[, 1], xx[, 2], z, phi=45, theta=45,
          xlab = "B_lf1", ylab = "hmat", zlab = "log(MAE)")
  }
  
  pred <- predGPsep(gp_ptr, xx, lite=TRUE)
  
  plot(xx[, par], pred$mean, type = "l", ylim = c(0, max(evals$y)), 
       xlab = "B_lf1", ylab = "log(MAE)")

  lines(xx[, par], pred$mean + 1.96 * sqrt(pred$s2), col=2, lty=2)
  lines(xx[, par], pred$mean - 1.96 * sqrt(pred$s2), col=2, lty=2)
  
  points(evals$x[, par], evals$y, col = "red")
  
}
