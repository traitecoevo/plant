## Support code:
make_equilibrium_runner <- function(p, ctrl) {
  pretty_num_collapse <- function(x, collapse=", ") {
    paste0("{", paste(prettyNum(x), collapse=collapse), "}")
  }
  
  p <- validate(p)
  
  # default is about 10 ind.m-2
  large_offspring_arriving_change <- ctrl$equilibrium_large_birth_rate_change
  
  # traverse list of strategies and pull birth_rates
  last_arrival_rate <- purrr::map_dbl(p$strategies, 
                                      ~ purrr::pluck(., birth_rate_y))
  
  default_schedule_times <- rep(list(p$node_schedule_times_default),
                                length(last_offspring_arriving))
  
  last_schedule_times <- p$node_schedule_times
  history <- NULL
  counter <- 1L
  
  function(birth_rate) {
    
    # if a runner has diverged significantly in the last iteration then
    # reset the schedule to it's defaults and rebuild
    if (any(abs(birth_rate - last_arrival_rate) > large_offspring_arriving_change)) {
      p$node_schedule_times <- default_schedule_times
    }
    
    # set birth rates
    p$strategies <- purrr::map(p$strategies, ~ .$birth_rate_y = birth_rate)
    
    # update schedule - starts from prev. schedule so should be fast for fine scale resolution
    p_new <- build_schedule(p, ctrl = ctrl)
    
    # TODO: change behaviour of `build_schedule` to objects rather than attributes
    offspring_production <- attr(p_new, "offspring_production", exact=TRUE)
    
    
    # (ANDREW) Double arrow modifies counter in parent environment.. 
    # I don't love this approach to counting, as it introduces a side effect
    # to an otherwise functional programming design 
    
    ## These all write up to the containing environment:
    p <<- p_new
    last_schedule_times <<- p_new$node_schedule_times
    history <<- c(history, list(c("in" = birth_rate, "out" = offspring_production)))
    
    msg <- sprintf("eq> %d: %s -> %s (delta = %s)", counter,
                   pretty_num_collapse(birth_rate),
                   pretty_num_collapse(offspring_production),
                   pretty_num_collapse(offspring_production - birth_rate))
    plant_log_eq(msg,
                 stage="runner",
                 iteration=counter,
                 birth_rate=birth_rate,
                 offspring_production=offspring_production)
    
    
    counter <<- counter + 1L
    
    # TODO: check why returning node schedule as an attribute
    # attr(offspring_production, "schedule_times") <- last_schedule_times
    
    offspring_production
  }
}

equilibrium_runner_cleanup <- function(runner, converged=TRUE) {
  ## This is super gross.
  e <- environment(runner)
  if (is.function(e$runner_full)) {
    runner <- e$runner_full
    e <- environment(runner)
  }
  
  p <- e$p
  
  # f <- function(s, br){
  #   s$birth_rate_y <- br
  #   return(s)
  # }
  
  #this may or may not work (original function below)     
  # p$birth_rate <- as.numeric(e$last_offspring_arriving)
  
  # browser()
  # p$strategies <- mapply(f, p$strategies, as.numeric(e$last_offspring_arriving), SIMPLIFY = FALSE)
  
  p$node_schedule_times <- e$last_schedule_times
  attr(p, "progress") <- rbind_list(e$history)
  attr(p, "converged") <- converged
  p
}

##' Check low-abundance strategies for viability.
##'
##' @title Check low-abundance strategies for viability
##' @param p A Parameters object
##' @param ctrl Control object
##' @export
check_inviable <- function(p, ctrl) {
  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a offspring arrival of less than
  ## `eps_test * max(p$birth_rate)` will be tested.  By default
  ##  this is 1 100th of the maximum offspring arrival.
  ## TODO: don't do anything if we don't have at least 2 species?
  eps <- ctrl$equilibrium_extinct_offspring_arriving
  ## TODO: This was ctrl$equilibrium_inviable_test, but I think
  ## that birth offspring arrival actually makes more sense?  It's fractional
  ## though so who knows.
  eps_test <- 1e-2
  birth_rate <- sapply(p$strategies, function(s) s$birth_rate_y, simplify = TRUE)
  ## NOTE: We don't actually run to equilibrium here; this is just
  ## because it's a useful way of doing incoming -> outgoing offspring
  ## rain.
  runner <- make_equilibrium_runner(p,ctrl =ctrl)
  offspring_production <- runner(offspring_production)
  
  test <- which(offspring_production < birth_rate &
                  birth_rate < max(offspring_production) * eps_test)
  test <- test[order(offspring_production[test])]
  
  drop <- logical(length(offspring_production))
  
  for (i in test) {
    plant_log_inviable(paste("Testing species", i),
                       stage="testing", species=i)
    x <- offspring_production
    x[i] <- eps
    res <- runner(x)
    if (res[[i]] < eps) {
      plant_log_inviable(paste("Removing species", i),
                         stage="removing", species=i)
      drop[[i]] <- TRUE
      res[[i]] <- 0.0
      offspring_production <- res
    }
  }
  
  ## It's possible that things slip through and get driven extinct by
  ## the time that they reach here.
  drop <- drop | offspring_production < eps
  
  attr(offspring_production, "drop") <- drop
  offspring_production
}