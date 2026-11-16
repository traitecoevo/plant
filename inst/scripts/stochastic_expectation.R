## Is the deterministic solver the large-patch limit of the stochastic one?
##
## With competition the two solvers need not agree on a finite patch: the
## environment is a nonlinear function of a random population, so a small
## patch's average stand differs from the deterministic one. As the patch grows
## the population's relative noise shrinks and the gap should close. This
## script measures the gap at several patch ages and areas, for each model's
## default strategy, and reports it in the stochastic mean's standard errors.
## A gap that stays put as the area grows is a disagreement between the
## solvers, not noise. (Without competition they must agree exactly, which
## tests/testthat/test-stochastic-expectation.R checks.)
##
## Run against the installed plant:
##
##   Rscript inst/scripts/stochastic_expectation.R [results.csv] [quick]
##
## It writes one row per model, patch age, area and death-check interval:
##
##   quantity     plants per m2 (density), and per m2 totals of height and of
##                competition_effect (what each plant contributes to the
##                environment: leaf area in FF16 and TF24, k_I times basal
##                area in K93);
##                for TF24 also the mean soil moisture over its layers
##   det          the deterministic value, cohorts every 0.05 yr
##   det_error    how far it moves from cohorts twice as far apart, a bound on
##                its own error where it is converging. Where det_error is
##                not small against det, the deterministic answer has not
##                converged and the gap says nothing about the stochastic one.
##   mean, se     the stochastic mean over `reps` independent patches
##   rel_gap      (mean - det) / det
##   z            (mean - det) / se; with as few as 20 patches it follows a t
##                distribution on reps - 1 degrees of freedom, not a normal
##
## Each stochastic patch is its own run ending at the patch age, which is a
## death-check time. Replicates are sized so that every area holds about the
## same number of plants in total, and run in parallel, each with its own seed.
## "quick" runs a small version, to check the script runs.
##
## One known difference between the solvers is not a matter of patch size. The
## stochastic patch sums the competition above a height over its plants oldest
## first and stops at the first that is not taller, which is exact only while
## older plants are taller. Where a model's cohorts cross, as TF24's can, a gap
## that does not close with area may be that approximation.

suppressPackageStartupMessages(library(plant))

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1L) args[[1L]] else "stochastic_expectation.csv"
quick <- "quick" %in% args
cores <- max(1L, parallel::detectCores() - 1L)

settings <- list(
  FF16 = list(ages = c(5, 10, 20, 40), areas = c(5, 20, 80, 320),
              plants = 4000, intervals = c(0.05, 0.01)),
  K93 = list(ages = c(5, 10, 20, 40), areas = c(5, 20, 80, 320),
             plants = 4000, intervals = c(0.05, 0.01)),
  TF24 = list(ages = c(5, 10), areas = c(5, 20, 80),
              plants = 2000, intervals = 0.05))
if (quick) {
  settings <- lapply(settings, function(s) {
    s$ages <- head(s$ages, 2)
    s$areas <- head(s$areas, 2)
    s$plants <- 200
    s$intervals <- s$intervals[1]
    s
  })
}

## TF24's soil cannot be integrated stochastically at Control()'s largest step;
## the test suite caps it the same way (see stochastic_control() there).
model_control <- function(x) {
  ctrl <- Control()
  if (x == "TF24") {
    ctrl$ode_step_size_max <- 0.05
  }
  ctrl
}

model_parameters <- function(x, t_max) {
  p <- base_parameters(x)
  p$strategies <- list(get(paste0(x, "_Model"))())
  p$max_patch_lifetime <- t_max
  p
}

## The stand of a deterministic run at each of `ages`, per m2.
deterministic_stand <- function(x, ages, spacing) {
  t_max <- max(ages)
  p <- model_parameters(x, t_max)
  ## Integer multiples, so each age is a schedule time exactly. The last
  ## cohort is one interval short of the end, which has no step of its own.
  p$node_schedule_times <- list(seq_len(round(t_max / spacing)) * spacing - spacing)
  ## Carried in birth date, the stand is read off the cohorts by integrating
  ## over birth dates, which never cross. Over height it is not to be trusted
  ## where cohorts overtake one another, as TF24's do.
  ctrl <- model_control(x)
  ctrl$node_density_coordinate <- "birth_date"
  res <- run(p, control = ctrl, record = "steps")
  td <- tidy(res)
  stand <- integrate_over_size_distribution(td$species)
  out <- lapply(ages, function(a) {
    row <- stand[abs(stand$time - a) < 1e-9, ]
    if (nrow(row) != 1L) {
      stop("No single step at age ", a)
    }
    v <- c(density = row$density, height = row$height,
           competition_effect = row$competition_effect)
    if (!is.null(td$env$soil_moist)) {
      v[["soil_moist"]] <- mean_soil_moist_at(td$env$soil_moist, a)
    }
    v
  })
  names(out) <- ages
  out
}

## The mean over soil layers at time `a`, from the last step there.
mean_soil_moist_at <- function(soil, a) {
  at <- soil[abs(soil$time - a) < 1e-9, ]
  mean(at$soil_moist[at$step == max(at$step)])
}

## One stochastic patch of `area` m2, run to `age`: its stand per m2.
stochastic_stand <- function(x, age, area, interval, seed) {
  p <- model_parameters(x, age)
  p$patch_area <- area
  ctrl <- model_control(x)
  ctrl$stochastic_death_interval <- interval
  set.seed(seed)
  res <- run(p, control = ctrl, method = "stochastic")
  st <- res$solver$patch$state
  plants <- st$species[[1]]
  alive <- attr(plants, "is_alive")
  v <- c(density = sum(alive),
         height = sum(plants["height", alive]),
         competition_effect = sum(plants["competition_effect", alive])) / area
  if (!is.null(st$env$soil_moist)) {
    v[["soil_moist"]] <- mean(st$env$soil_moist)
  }
  v
}

rows <- list()
seed <- 0L
for (x in names(settings)) {
  s <- settings[[x]]
  message(x, ": deterministic")
  det <- deterministic_stand(x, s$ages, 0.05)
  det_coarse <- deterministic_stand(x, s$ages, 0.1)
  for (age in s$ages) {
    d <- det[[as.character(age)]]
    d_error <- abs(det_coarse[[as.character(age)]] - d)
    for (area in s$areas) {
      reps <- min(500L, max(20L, ceiling(s$plants / (area * d[["density"]]))))
      for (interval in s$intervals) {
        message(sprintf("%s: age %g, area %g, interval %g, %d patches",
                        x, age, area, interval, reps))
        seeds <- seed + seq_len(reps)
        seed <- seed + reps
        runs <- parallel::mclapply(seeds, function(sd)
          stochastic_stand(x, age, area, interval, sd), mc.cores = cores)
        ## A run that threw comes back as a try-error, and one whose worker
        ## was killed (out of memory, say) as NULL, which rbind() would drop
        ## while the standard error still counted it.
        failed <- vapply(runs, function(r) is.null(r) || inherits(r, "try-error"),
                         logical(1))
        if (any(failed)) {
          first <- runs[[which(failed)[1]]]
          stop(sum(failed), " stochastic runs failed: ",
               if (is.null(first)) "a worker was lost"
               else conditionMessage(attr(first, "condition")))
        }
        m <- do.call(rbind, runs)
        for (q in names(d)) {
          mean_q <- mean(m[, q])
          se_q <- stats::sd(m[, q]) / sqrt(reps)
          rows[[length(rows) + 1L]] <- data.frame(
            model = x, age = age, area = area, interval = interval,
            reps = reps, quantity = q, det = d[[q]], det_error = d_error[[q]],
            mean = mean_q, se = se_q, rel_gap = (mean_q - d[[q]]) / d[[q]],
            z = (mean_q - d[[q]]) / se_q)
        }
      }
    }
  }
}

results <- do.call(rbind, rows)
utils::write.csv(results, out_file, row.names = FALSE)
print(within(results, {
  det <- signif(det, 4); det_error <- signif(det_error, 2)
  mean <- signif(mean, 4); se <- signif(se, 2)
  rel_gap <- signif(rel_gap, 2); z <- round(z, 1)
}), row.names = FALSE)
message("Written to ", out_file)
