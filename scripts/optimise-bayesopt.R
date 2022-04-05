library(zeallot)
library(glue)
library(GPfit)
library(animation)

# Source: https://github.com/bearloga/bayesopt-tutorial-r/blob/master/BayesOpt.R

# Goal: find x and y which maximize f
f <- function(x) {
  return((6 * x - 2)^2 * sin(12 * x - 4))
}

# generate starting points
c(max_evals, seed_evals) %<-% c(8, 4)

# matrix to store evaluations of f:
evaluations <- matrix(
  as.numeric(NA),
  ncol = 2, nrow = seed_evals,
  dimnames = list(NULL, c("x", "y"))
)

# seed with a few evaluations:
evaluations[1:seed_evals, "x"] <- seq(0, 1, length.out = seed_evals)
evaluations[1:seed_evals, "y"] <- f(evaluations[1:seed_evals, "x"])

# plot f(x) against f
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1))
curve(f(x), 0, 1)
points(evaluations[1:seed_evals, ], pch = 16)


# BayesOpt tells us where to sample next
set.seed(0)
bayesian_optimize <- function(optmize_function, 
                              init_evals, 
                              max_iter, 
                              acquisition_function,
                              minimize = TRUE, 
                              control = NULL) {
  
  # expand to hold additional evaluations:
  evaluations <- rbind(init_evals, matrix(
    as.numeric(NA),
    ncol = 2, nrow = max_iter,
    dimnames = list(NULL, c("x", "y"))
  ))
  
  if (is.null(control)) {
    control <- list(cov = list(type = "exponential", power = 1.95))
    # control <- list(type = "matern", nu = 5/2)
  }
  if (acquisition_function == "cb") {
    if (is.null(control$kappa)) {
      control$kappa <- 2
    }
  }
  x_new <- seq(0, 1, length.out = 100) # potential x's to evaluate
  
  for (eval_iter in (nrow(init_evals) + 1):(nrow(init_evals) + max_iter)) {
    
    fit <- GP_fit(
      X = evaluations[1:(eval_iter - 1), "x"],
      Y = evaluations[1:(eval_iter - 1), "y"],
      corr = control$cov
    )
    
    predictions <- predict.GP(fit, xnew = data.frame(x = x_new))
    mu <- predictions$Y_hat
    sigma <- sqrt(predictions$MSE)
    
    if (minimize) {
      y_best <- min(evaluations[, "y"], na.rm = TRUE)
    } else {
      y_best <- max(evaluations[, "y"], na.rm = TRUE)
    }
    
    
    # three aq. functions: 
    # poi - Probability of improvement
    # ei - expected improvement
    # cb - upper/lower confidence bound
    
    if (acquisition_function == "poi") {
    
      # Probability of improvement:
      acquisition <- purrr::map2_dbl(mu, sigma, function(m, s) {
        if (s == 0) return(0)
        else return(pnorm((y_best - m) / s))
      })
      if (!minimize) {
        acquisition <- 1 - acquisition
      }
      
      x_next <- x_new[which.max(acquisition)]
      plot(x_new, acquisition, type = "l", col = "red", 
           ylim = c(0, 1), xlab = "x", 
           ylab = expression("a"["POI"]))
      
    } else if (acquisition_function == "ei") {
      
      # Expected improvement:
      acquisition <- purrr::map2_dbl(mu, sigma, function(m, s) {
        if (s == 0) return(0)
        gamma <- (y_best - m) / s
        if (minimize) {
          phi <- pnorm(gamma)
        } else {
          phi <- 1 - pnorm(gamma)
        }
        return(s * (gamma * phi + dnorm(gamma)))
      })
      
      x_next <- x_new[which.max(acquisition)]
      plot(x_new, acquisition, type = "l", col = "red",
           xlab = "x", ylab = expression("a"["EI"]))
      
    } else if (acquisition_function == "cb") {
      
      # GB upper/lower confidence bound:
      if (minimize) {
        acquisition <- mu - control$kappa * sigma
        x_next <- x_new[which.min(acquisition)]
        plot(x_new, acquisition, type = "l", col = "red", 
             xlab = "x", ylab = expression("a"["LCB"]))
        
      } else {
        acquisition <- mu + control$kappa * sigma
        x_next <- x_new[which.max(acquisition)]
        plot(x_new, acquisition, type = "l", col = "red", 
             xlab = "x", ylab = expression("a"["UCB"]))
      }
    } else {
      stop("acquisition_function must be 'poi', 'ei', 'cb'")
    }
    
    abline(v = x_next, lty = "dashed", col = "red", lwd = 2)
    acquisition_function_label <- switch(
      acquisition_function,
      "poi" = "probability of improvement",
      "ei" = "expected improvement",
      "cb" = paste("GP", ifelse(minimize, "lower"), "confidence bound")
    )
    
    legend("topleft", glue("proposal via {acquisition_function_label}"),
           bty = "n", col = "red", lty = "dashed", lwd = 2)
    
    # Visualize hidden function and GP fit:
    curve(f(x), 0, 1, lwd = 1.5)
    lines(x_new, mu, col = "blue", lwd = 2, lty = "dotted")
    polygon(c(x_new, rev(x_new)), c(mu + sigma, rev(mu - sigma)),
            col = rgb(0, 0, 1, 0.25), border = NA)
    
    points(evaluations, pch = 16)
    points(evaluations[(eval_iter - 1), "x"], 
           evaluations[(eval_iter - 1), "y"],
           pch = 16, col = "red")
    
    abline(v = x_next, lty = "dashed", col = "red", lwd = 2)
    legend("topleft", "most recent evaluation", bty = "n", col = "red", pch = 16)
    
    y_next <- f(x_next)
    evaluations[eval_iter, ] <- c(x_next, y_next)
  }
  
  return(list(x = x_next, y = y_next))
}

bayesian_optimize(f, evaluations, max_evals - seed_evals, "ei")


# Visualize | mar = c(5.1, 4.1, 4.1, 2.1) (bottom, left, top, right)
purrr::walk(c("poi", "ei", "cb"), function(af) {
  # Static images:
  png(glue("bayesopt_{af}.png"), width = 12, height = 12, units = "in", res = 300)
  par(mfrow = c(4, 2), mar = c(4.1, 4.1, 0.5, 0.5), cex = 1.1)
  bayesian_optimize(f, evaluations, max_evals - seed_evals, af)
  dev.off()
  # Animated GIF:
  saveGIF(
    {
      par(mfrow = c(1, 2), mar = c(4.1, 4.1, 0.5, 0.5), cex = 1.1)
      bayesian_optimize(f, evaluations, max_evals - seed_evals, af)
    },
    glue("bayesopt_{af}.gif"), nmax = 4, loop = TRUE, interval = 1.5,
    ani.width = 900, ani.height = 300, ani.dev = "png",
    autobrowse = FALSE
  )
})
