library(tree)
library(testthat)

p <- new(Parameters)
p$set_parameters(list(patch_area=1.0))   # See issue #13
p$set_control_parameters(fast_control()) # A bit faster

hmat_r <- c(1, 50) # only fails when lower bound is 1, not 2.
trait <- "hmat"
values <- seq_log(hmat_r[[1]], hmat_r[[2]], 20)

## Set up a schedule based on the middle strategy:
values_mid <- exp(mean(log(range(values))))
p1 <- expand_parameters(trait, values_mid, p, mutant=FALSE)
ebt <- run_ebt(p1, default_cohort_schedule(p1))
schedule <- ebt$cohort_schedule

## Then just add a single problem-causing strategy -- this is just a
## strategy that has a height much shorter than the one used to
## establish the ode times.
value <- 1
p_empty <- p$copy()
p_empty$clear()
p_empty <- expand_parameters(trait, value, p_empty, mutant=TRUE)

schedule_empty <- new(CohortSchedule, 1)
schedule_empty$max_time  <- schedule$max_time
schedule_empty$ode_times <- schedule$ode_times
schedule_empty$all_times <-
  rep(list(unique(sort(unlist(schedule$all_times)))), 1)
schedule_empty$use_ode_times <- TRUE

## This generates an error.
expect_that(ebt_empty <- run_ebt(p_empty, schedule_empty),
            throws_error("height must be positive"))

## Run the system up until the error:
ebt <- new(EBT, p_empty)
ebt$cohort_schedule <- schedule_empty
try(ebt$run())

## Extract the patch at this point:
patch <- ebt$patch

m <- matrix(ebt$ode_values, 4)
rownames(m) <- c("height", "mortality", "fecundity", "log_density")
r <- matrix(ebt$ode_rates, 4)
rownames(r) <- c("height", "mortality", "fecundity", "log_density")

## Here are the cohort heights, showing the instability kick in for
## later cohorts as the step size increases (this should be monotonic)
plot(m["height",])
