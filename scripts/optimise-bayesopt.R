# source: https://bookdown.org/rbg/surrogates/chap7.html
library(laGP)
library(lhs)

eps <- sqrt(.Machine$double.eps) ## used lots below


# strategy defaults
mulga <- function() {
  p0 <- scm_base_parameters("FF16bg", "FF16_Env")
  
  p0$strategy_default$lma <- 0.0645
  p0$strategy_default$hmat <- 5
  p0$strategy_default$rho <- 1100
  p0$strategy_default$narea <- 0.0032
  p0$strategy_default$a_l1 <- 2.17
  p0$strategy_default$a_l2 <- 0.5
  p0$strategy_default$omega <- 0.00000926
  p0$strategy_default$k_2 <- 0.00
  p0
}


create_environment <- function() FF16bg_make_environment()
create_site <- make_FF16bg_hyperpar

create_species <- function(lma = 0.0645, hmat = 5, eta = 12, k_2 = 0,
                           birth_rate = 100,
                           site) {
  
  traits = trait_matrix(c(lma, hmat, eta, k_2), c("lma", "hmat", "eta", "k_2"))
  sp = expand_parameters(traits, mulga(), site, mutant = FALSE,
                        birth_rate_list = list(birth_rate))

  # base <- scm_base_parameters("FF16bg", "FF16_Env")
  # sp = expand_parameters(traits, base, site, mutant = FALSE, 
  #                      birth_rate_list = list(birth_rate))
  # 
   
  return(sp)
}

# This is clunky - assumes one schedule will be suitable even while optimising
# traits and site characteristics
create_schedule <- function(max_patch_lifetime = 250,
                            optimise_schedule = FALSE,
                            schedule_reduction_factor = 5) {
  
  site <- create_site()
  parameters <- create_species(site = site)
  
  # update patch longevity
  parameters$max_patch_lifetime <- max_patch_lifetime
  
  if(optimise_schedule) {
    
    # start with built-in schedule
    nodes <- node_schedule_times_default(max_patch_lifetime)
    
    # then downsample, taking every nth integration node 
    nth_element <- function(vector, n = 1, starting_position = 1) { 
      vector[c(1:starting_position, seq(starting_position, length(vector), n))] 
    }
    
    times <- nth_element(nodes, schedule_reduction_factor)
    
    parameters$node_schedule_times[[1]] <- times #1spp
    
    # and rebuild
    parameters <- build_schedule(parameters)
  }
  
  return(parameters)
}


calculate_net_biomass <- function(parameters, environment = create_environment()) {
  
  # gather outputs at each time step
  results <- run_scm_collect(parameters, environment) %>% 
    tidy_patch() %>%
    FF16_expand_state()
  
  v <- c("mass_leaf", "mass_bark", "mass_sapwood", "mass_heartwood")
  
  tissues <- results$species %>% 
    integrate_over_size_distribution() %>%
    select(time,species, one_of(v)) %>%
    pivot_longer(cols=starts_with("mass"), names_to = "tissue") %>%
    mutate(across(tissue, factor, levels = v)) %>%
    filter(time != max(time)) # drop last node
  
  return(tissues)
}



# @param x optimisable parameters
# @param target a function or vector of points to calculate the marginal likelihood 
# @param latitude used to specify site characteristics
# @param n_control_points number of points to evaluate LL if target is a function
calibrate <- function(x, target, 
                      latitude = 28.182,
                      parameters, 
                      environment = NULL,
                      match_tissue = "total_biomass",
                      plot_comparison = TRUE) {
  
  if(is.null(parameters$node_schedule_times)) {
    stop("Please run `create_schedule` first")
  }
  
  # assume optimising B_lf1 and hmat
  with(as.list(x), {
    site <- create_site(B_lf1 = B_lf1, 
                        latitude = latitude)
    
    sp <- create_species(hmat = hmat, 
                         eta = eta, 
                         k_2 = k_2,
                         site = site)

    # overwrite species
    parameters$strategies <- sp$strategies

    if(is.null(environment)) {
      environment = create_environment()
    }
    
    biomass = calculate_net_biomass(parameters, environment) 
    
    if(is.function(target)) {

      if(match_tissue == "total_biomass") {
        biomass <- group_by(biomass, time, species) %>%
          summarise(tissue = "total_biomass",
                    value = sum(value), 
                    .groups = "drop")
      }
      
      # something odd with the last point
      predictions <- filter(biomass, tissue == match_tissue) 
      
      if(plot_comparison == T) {
        ggplot(predictions, aes(time, value)) +
          geom_area(aes(fill = tissue)) +
          geom_function(fun=tyf, linetype = "dashed") +
          labs(x = "Patch age (yr)", y = "Above ground mass (kg/m2)") +
          theme_classic() + 
          facet_wrap(~species)
      }
      
      # https://bookdown.org/rbg/surrogates/chap8.html#chap8calibopt
      # fit GP to plant results
      x <- matrix(predictions$time, ncol = 1)
      
      da <- darg(list(mle=TRUE), x)
      ga <- garg(list(mle=TRUE), predictions$value)
      
      surrogate <- newGPsep(x, predictions$value, d=da$start, g=ga$start, dK=TRUE)
      mleGPsep(surrogate, 
               param = "d", 
               tmin=c(da$min),
               tmax=c(da$max))
      
      z <- predGPsep(surrogate, x, lite=T, nonug=T)$mean
      

      # fit GP to residuals
      delta <- target(x) - z
      residuals <- newGPsep(x, delta, d=da$start, g=ga$start, dK=TRUE)
      
      mleGPsep(residuals, 
               param = "both", 
               tmin=c(da$min, ga$min),
               tmax=c(da$max, ga$max),
               ab = c(da$ab, ga$ab))

      # plot to see interpolation of residuals 
      if(FALSE) {
        xx = matrix(seq(0, max(x), len = 100), ncol=1)
        plot(xx, predGPsep(residuals, xx, lite=T, nonug=T)$mean)
        points(x, delta, col = "red", type = "l")
      }
      
      # extract neg. log. lik.
      nll <- -llikGPsep(residuals, dab = da$ab, gab = ga$ab)
      
      deleteGPsep(surrogate)
      deleteGPsep(residuals)
      
    } else {
      stop("Calibration for dataframes not yet implemented")
    }
    
    return(nll)
  })
}


# https://bookdown.org/rbg/surrogates/chap7.html#classic-ei-illustration
EI.search <- function(x, y, gp, bounds, multi_start = 5, tol = eps) {

  # find minima of data and start there
  m <- which.min(y)
  fmin <- y[m]
  start <- matrix(x[m, ], nrow=1)
  
  # optionally add extra samples
  if(multi_start > 1) {
    n = multi_start - 1
    jitter <- generate_seeds(gp, bounds, n = n, evaluate = F)$x_raw
    start <- rbind(start, jitter)
  } else {
    message("Using a multi-start optimisation routine can improve numerical stability")
  }
  
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
  
  max_ei = 0

  # optim likes to minimise things, in this case the negative expected 
  # improvement of the next acquisition of the neg. log. likelihood (curly!)
  objective <- function(x, fmin, gp) {
    - EI(gp, x, fmin)
  }
  
  # repeatedly optimise
  for(i in 1:nrow(start)) {
    
    if(EI(gp, start[i, ], fmin) <= tol) {
      out <- list(value = -Inf)
      next 
    }

    out <- optim(start[i, ], objective, method = "L-BFGS-B", 
                 lower = rep(0, length(bounds)),
                 upper = rep(1, length(bounds)),
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


# https://bookdown.org/rbg/surrogates/chap7.html#chap7ieci
IECI.search <- function(x, y, gp, bounds, multi_start = 5, tol = eps) {
  
  # find minima of data and start there
  fmin <- min(predGPsep(gp, x, lite=TRUE, nonug = T)$mean)

  m <- which(y == min(y))
  start <- x[m, ]

  # optionally add extra samples
  if(multi_start > 1) {
    n = multi_start - 1
    jitter <- generate_seeds(gp, bounds, n = n, evaluate = F)$x_raw
    start <- rbind(start, jitter)
  } else {
    message("Using a multi-start optimisation routine can improve numerical stability")
  }
  
  max_ieci = 0

  # optim likes to minimise things, in this case the negative expected 
  # conditional improvement of the next acquisition of the neg. log. likelihood (curly!)
  objective <- function(x_next, fmin, gp) {
    xx <- seq(0, 1, len = 100)
  
     -ieciGPsep(gp, matrix(x_next, ncol = length(bounds)), fmin, 
                Xref = matrix(rep(xx, length(bounds)),
                              ncol=length(bounds)), nonug=TRUE)
  }
  
  # repeatedly optimise
  for(i in 1:nrow(start)) {
    out <- optim(start[i, ], objective, 
                 method = "L-BFGS-B", 
                 lower = rep(0, length(bounds)),
                 upper = rep(1, length(bounds)),
                 gp = gp, fmin = fmin)
    
    ieci = -out$value
    
    # save best acquisition location
    if(ieci > max_ieci) {
      max_ieci = ieci
      x_next = matrix(out$par, nrow = 1)
    } 
  }
  
  return(x_next)
}



bayesopt <- function(f, target, bounds, parameters, n_iter = 1, n_init = 6,
                     gp = NULL, evals = NULL, search_fn = "EI") {

  if(is.null(bounds))
    stop("Please specify the bounds of optimisable parameters")
    
  if(is.null(evals)) {
    if(n_init < 1)
      stop("Provide a set of initial simulations (seed) or 
           randomly generate n_init new simulations")
    
    evals <- generate_seeds(f, n_init, target, bounds, parameters)
  }
  
  x_raw <- evals$x_raw
  y <- evals$y
  
  # initialise GP in C, with separable correlation for >1 dimensions
  # d = length scale, g = nugget, dK = save derivative for optimisation of d
  if(is.null(gp)) {
    gp_ptr <- newGPsep(x_raw, y, d=0.1, g=1e-6, dK=TRUE)
    
    # priors for GP lengthscale, optimised using MLE
    da <- darg(list(mle = T), 
               generate_seeds(calibrate, bounds, n = 1000, evaluate = F)$x_raw)

  } else if(is.list(gp)) {
    gp_ptr = gp$pointer
    da <- gp$priors
  } else {
    stop("GP should have a pointer and a priors object")
  }
  
  mleGPsep(gp_ptr, param="d", tmin= da$min, tmax= da$max)  
  
  ## optimization loop of sequential acquisitions
  for(i in 1:n_iter) {
    
    # find best expected conditional improvement (IECI)
    if(search_fn == "IECI") {
      search = IECI.search
    } else {
      search = EI.search
    }
    
    x_next_raw <- search(x_raw, y, gp_ptr, bounds)
    colnames(x_next_raw) <- names(bounds)

    # collect data at proposed location
    x_next <- unscale(x_next_raw, bounds)
    y_next <- apply(x_next, 1, f, target = target, parameters = parameters)
    
    # update GP and length-scale
    updateGPsep(gp_ptr, x_next_raw, y_next)
    mleGPsep(gp_ptr, param = "d", tmin = da$min, tmax = da$max, ab = da$ab)
    
    # save new simulation
    x_raw <- rbind(x_raw, x_next_raw)
    y <- c(y, y_next)
  }
  

  return(list(gp = list(pointer = gp_ptr, priors = da), 
              evaluations = list(x_raw = x_raw, y = y),
              minima = min(y)))
}


unscale <- function(x_raw, bounds) {
  x = x_raw
  
  for(i in 1:length(bounds)) {
    diff = bounds[[i]]$max - bounds[[i]]$min
    x[, i] = x_raw[, i] * diff + bounds[[i]]$min
  }
  
  return(x)
}

# @param bounds a list with upper and lower limits for each variable
generate_seeds <- function(f, bounds, n = 1,
                           target = NULL,
                           parameters = NULL,
                           evaluate = T) {
  x_raw <- randomLHS(n, length(bounds))
  colnames(x_raw) <- names(bounds)
  
  if(evaluate == T) {
    if(is.null(target) | is.null(parameters)) {
      stop("Target and parameters must be provided to evaluate seeds")
    }
    
    x <- unscale(x_raw, bounds)
    y <- apply(x, 1, f, target = target, parameters = parameters)
    return(list(x_raw = x_raw, y = y))
    
  } else {
    return(list(x_raw = x_raw))
  }
}


generate_seeds_parallel <- function(f, bounds, n = 1,
                           target = NULL,
                           parameters = NULL,
                           evaluate = T) {
  x_raw <- randomLHS(n, length(bounds))
  colnames(x_raw) <- names(bounds)
  
  if(evaluate == T) {
    if(is.null(target) | is.null(parameters)) {
      stop("Target and parameters must be provided to evaluate seeds")
    }
    
    x <- unscale(x_raw, bounds)
    
    y <- foreach(i=1:n, .combine=c) %dopar%
      f(x[i, ], target = target, parameters = parameters)
    
    return(list(x_raw = x_raw, y = y))
  } else {
    return(list(x_raw = x_raw))
  }
}



generate_range <- function(bounds, n = 1000) {
  x = matrix(nrow = n, ncol= length(bounds))
  basis = seq(0, 1, len = n)
  
  for(i in 1:length(bounds)) {
    x[, i] = basis
  }
  
  return(x)
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
  
  x <- unscale(xx, bounds)[, par]
  plot(x, pred$mean, type = "l", 
       xlab = "B_lf1", ylab = "negLogLik")

  lines(x, pred$mean + 1.96 * sqrt(pred$s2), col=2, lty=2)
  lines(x, pred$mean - 1.96 * sqrt(pred$s2), col=2, lty=2)
  
  
  points(unscale(evals$x, bounds)[, par], evals$y, col = "red")
  
}
