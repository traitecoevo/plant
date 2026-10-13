## A continuous patch run to its fixed point (#519; plan decision 9): constant
## drivers, no disturbance, a cohort every PLANT_519_DT years and the absorbing
## top layer. Four tasks, run in this order by default (PLANT_519_TASK picks
## some):
##
##   settle      block means until the stand settles, then on to
##               PLANT_519_VERIFY to say whether it is a fixed point or a
##               cycle (output/<model>-continuous.csv, -continuous-snapshots.csv)
##   cost        PLANT_519_TMAX years with the layer and without it; without
##               it every cohort stays integrated, so on TF24 this is slow
##   boundaries  the stand at three absorbing heights against the unbounded
##               one: cohorts integrated, the top of the stand, its census and
##               fitness (output/<model>-boundaries.csv)
##   landscape   the resident's equilibrium birth rate b*, and a mutant lma
##               sweep there, each mutant one cohort held in the resident's
##               field (fixed_point_fitness()), with three checks: the field
##               held at four points inside one interval, held a century
##               either side, and against the same cohort born into the
##               running stand (output/<model>-equilibrium-mutants.csv, -check.csv)
##
## ⚠️ THE INTRODUCTION GRID HAS TO BE FINE ENOUGH FOR A FIXED POINT. Between two
## introductions the field swings, so a grid that is too coarse makes the stand
## a periodic orbit, and a field held at an introduction is the top of the
## swing. For FF16 at b* the swing in held R is 1.8 % at a cohort a year, 0.17 %
## at two and 0.02 % at four (plan hazard 12); `landscape` uses four.
##
##   PLANT_519_MODEL     "FF16" (default) or "TF24"
##   PLANT_519_ABSORB    control$absorb_epsilon (0.02 for FF16, 0.05 for TF24)
##   PLANT_519_RAIN      TF24's constant rainfall, m per year (1)
##   PLANT_519_DT        years between introductions (1 for settle, cost and
##                       boundaries, 0.25 for landscape)
##   PLANT_519_BIRTH     birth rate for settle, cost and boundaries (1)
##   PLANT_519_TMAX      cap on settle, and the span of cost (3000)
##   PLANT_519_BLOCK     years per convergence block (50)
##   PLANT_519_TOL       relative change counted as settled (0.005)
##   PLANT_519_VERIFY    year settle runs on to after converging (1000)
##   PLANT_519_T_EQ      year at which landscape reads the field (400)
##   PLANT_519_LIFETIME  years each mutant runs; it must outlive its cohort (3000)
##   PLANT_519_LIFE      years the in-stand reference carries its cohort (400)
##   PLANT_519_LMA       landscape's mutant lma values
##   PLANT_519_OUT       output directory

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(here, "common.R"))

model    <- Sys.getenv("PLANT_519_MODEL", "FF16")
tf24     <- identical(model, "TF24")
absorb   <- num_env("PLANT_519_ABSORB", if (tf24) 0.05 else 0.02)
rain     <- num_env("PLANT_519_RAIN", 1)
birth    <- num_env("PLANT_519_BIRTH", 1)
t_max    <- num_env("PLANT_519_TMAX", 3000)
block    <- num_env("PLANT_519_BLOCK", 50)
tol      <- num_env("PLANT_519_TOL", 0.005)
t_verify <- num_env("PLANT_519_VERIFY", 1000)
t_eq     <- num_env("PLANT_519_T_EQ", 400)
lifetime <- num_env("PLANT_519_LIFETIME", 3000)
life     <- num_env("PLANT_519_LIFE", 400)
lma_grid <- vec_env("PLANT_519_LMA", c(0.04, 0.05, 0.06, 0.07, 0.0825, 0.1, 0.125, 0.15, 0.2))
tasks <- tasks_to_run(c("settle", "cost", "boundaries", "landscape"))
prefix <- tolower(model)
resident_lma <- 0.0825
stamp("plant %s; %s, resident lma %g, absorb_epsilon %g%s; tasks %s",
      as.character(packageVersion("plant")), model, resident_lma, absorb,
      if (tf24) sprintf(", rainfall %g m/yr", rain) else "", paste(tasks, collapse = ", "))

env_type <- if (tf24) "TF24_Env" else "FF16_Env"
hyperpar <- if (tf24) TF24_hyperpar else FF16_hyperpar
## TF24's soil ODE is stiff for an explicit stepper without a cap on the step.
make_control <- function(absorb_epsilon = absorb) {
  if (tf24) {
    Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb_epsilon,
            ode_step_size_max = 0.05)
  } else {
    Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb_epsilon)
  }
}
make_env <- function() {
  env <- Environment(model)
  if (tf24) env$extrinsic_drivers_set_constant("rainfall", rain)
  env
}
add_lma <- function(p, lma, b) {
  add_strategies(p, trait_matrix(lma, "lma"), hyperpar = hyperpar, birth_rate = b)
}
## The resident, a cohort every `dt` years from time zero.
resident_parameters <- function(b, t_max, dt) {
  p <- scm_base_parameters(model)
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- t_max
  p <- add_lma(p, resident_lma, b)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, t_max, dt))
  p
}

## ---- settle ---------------------------------------------------------------------
if ("settle" %in% tasks) {
  dt <- num_env("PLANT_519_DT", 1)
  p <- resident_parameters(birth, max(t_max, t_verify), dt)
  scm <- run_scm(p, make_env(), make_control(), run = FALSE)
  scm$snapshot_interval <- 5
  rows <- list()
  previous <- NULL
  streak <- 0L
  converged <- FALSE
  k <- 0L
  repeat {
    k <- k + 1L
    t_end <- k * block
    if (t_end > t_max) break
    w <- system.time(scm$run_to(t_end))[["elapsed"]]
    sn <- tidy_snapshots(scm)$species
    in_block <- sn$time > t_end - block & sn$time <= t_end
    m <- c(density = mean(sn$density[in_block]), height_max = mean(sn$height_max[in_block]),
           n_live = mean(sn$n_live[in_block]))
    rel <- if (is.null(previous)) rep(NA_real_, 3) else abs(m - previous) / previous
    sp <- scm$patch$species[[1]]
    rows[[k]] <- data.frame(time = scm$time, wall_s = w, s_per_year = w / block,
                            steps = length(scm$ode_times), nodes = sp$size,
                            retired = sum(sp$retired), live = sp$size - sum(sp$retired), t(m),
                            rel_density = rel[1], rel_height_max = rel[2], rel_n_live = rel[3])
    if (k %% 4 == 0 || k <= 2) {
      stamp("%5.0f yr: %.3f s/yr, %d cohorts, %d live; block means density %.4g, height_max %.4g, live %.1f; rel change %s",
            scm$time, w / block, sp$size, sp$size - sum(sp$retired), m[1], m[2], m[3],
            paste(signif(rel, 2), collapse = "/"))
    }
    streak <- if (!is.null(previous) && isTRUE(all(rel < tol))) streak + 1L else 0L
    previous <- m
    if (streak >= 2L) {
      converged <- TRUE
      break
    }
  }
  d <- do.call(rbind, rows)
  stamp("%s at %g yr; wall %.1f s", if (converged) "CONVERGED" else "not converged", scm$time,
        sum(d$wall_s))
  ## Fixed point or cycle? On past convergence, and the snapshots from there.
  t_conv <- scm$time
  scm$run_to(t_verify)
  sn <- tidy_snapshots(scm)$species
  stat <- sn[sn$time >= t_conv, ]
  stamp("from %g to %g yr: density mean %.6g, CV %.2e; tallest live cohort %.4f to %.4f m; live cohorts %d to %d",
        t_conv, scm$time, mean(stat$density), sd(stat$density) / mean(stat$density),
        min(stat$height_max), max(stat$height_max), min(stat$n_live), max(stat$n_live))
  ## Where the density sits: by cohort age.
  sp <- scm$patch$species[[1]]
  age <- scm$time - sp$node_times
  nu <- exp(sp$log_densities_state)
  live <- nu > 0
  stamp("at %g yr: %d live cohorts, ages %g to %g; share of density in cohorts older than 50 / 100 yr: %.2e / %.2e; tallest live %.3f m",
        scm$time, sum(live), min(age[live]), max(age[live]), sum(nu[age > 50]) / sum(nu),
        sum(nu[age > 100]) / sum(nu), max(sp$heights[live]))
  write.csv(d, file.path(out_dir, paste0(prefix, "-continuous.csv")), row.names = FALSE)
  write.csv(sn, file.path(out_dir, paste0(prefix, "-continuous-snapshots.csv")), row.names = FALSE)
}

## ---- cost -----------------------------------------------------------------------
## What the top layer buys: the same span with and without it.
if ("cost" %in% tasks) {
  p <- resident_parameters(birth, t_max, num_env("PLANT_519_DT", 1))
  for (r in c(absorb, 0)) {
    s <- run_scm(p, make_env(), make_control(r), run = FALSE)
    w <- system.time(s$run())[["elapsed"]]
    sp <- s$patch$species[[1]]
    stamp("%g yr with absorb_epsilon %g: %.1f s, %d cohorts, %d still integrated, %d steps",
          t_max, r, w, sp$size, sp$size - sum(sp$retired), length(s$ode_times))
  }
}

## ---- boundaries ----------------------------------------------------------------
if ("boundaries" %in% tasks) {
  dt <- num_env("PLANT_519_DT", 1)
  rows <- list()
  for (eps in c(0, 0.05, 0.02, 0.01)) {
    s <- run_scm(resident_parameters(birth, t_verify, dt), make_env(), make_control(eps),
                 run = FALSE)
    s$snapshot_interval <- 50
    w <- system.time(s$run())[["elapsed"]]
    sp <- s$patch$species[[1]]
    late <- tidy_snapshots(s)$species
    late <- late[late$time >= t_verify / 2, ]
    census <- stand_census(s)
    rows[[length(rows) + 1]] <- data.frame(
      absorb_epsilon = eps, height = if (eps > 0) sp$top_layer$height else NA,
      wall_s = w, integrated = sp$size - sum(sp$retired),
      top_min = min(late$height_max), top_max = max(late$height_max),
      leaf_area = census[["leaf_area"]], biomass = census[["mass_above_ground"]],
      R = s$net_reproduction_ratios)
    stamp("absorb_epsilon %-5g (%s m): %.2f s; integrated %4d; from %g yr the top runs %.5f to %.5f m; leaf area %.9g, biomass %.9g, R %.9g",
          eps, if (eps > 0) format(sp$top_layer$height, digits = 6) else "none", w,
          sp$size - sum(sp$retired), t_verify / 2, min(late$height_max), max(late$height_max),
          census[["leaf_area"]], census[["mass_above_ground"]], s$net_reproduction_ratios)
  }
  d <- do.call(rbind, rows)
  ## Against the unbounded stand.
  for (col in c("leaf_area", "biomass", "R")) {
    d[[paste0("rel_", col)]] <- d[[col]] / d[[col]][d$absorb_epsilon == 0] - 1
  }
  print(format(d, digits = 6), row.names = FALSE)
  write.csv(d, file.path(out_dir, paste0(prefix, "-boundaries.csv")), row.names = FALSE)
}

## ---- landscape --------------------------------------------------------------
if ("landscape" %in% tasks) {
  dt <- num_env("PLANT_519_DT", 0.25)
  ## The resident run to `t_read` and stopped there, whether or not it is an
  ## introduction: the run ends at its lifetime, which is how a field is read
  ## inside an interval.
  resident_at <- function(b, t_read = t_eq) {
    p <- resident_parameters(b, t_read, dt)
    scm <- run_scm(p, make_env(), make_control(), run = FALSE)
    scm$run()
    list(scm = scm, p = p)
  }
  fitness <- function(r, lma) {
    fixed_point_fitness(r$scm, r$p, trait_matrix(lma, "lma"), hyperpar = hyperpar,
                        lifetime = lifetime)
  }
  ## The same cohort born into the running stand mid-interval, at a birth rate
  ## 1e-6 of the resident's; a window around it, holding no resident cohort,
  ## exempts it alone from absorption.
  in_stand <- function(b, lma) {
    p <- resident_parameters(b, t_eq + life, dt)
    p <- add_lma(p, lma, b * 1e-6)
    p$node_schedule_times[[2]] <- t_eq + 0.5 * dt
    p <- fitness_window(p, t_eq + 0.25 * dt, t_eq + 0.75 * dt)
    w <- system.time(scm <- run_scm(p, make_env(), make_control()))[["elapsed"]]
    sp <- scm$patch$species
    c(R_mutant = sp[[2]]$net_reproduction_ratio_by_node[[1]] * p$strategies[[2]]$pars$S_D,
      wall = w)
  }

  ## The resident's demographic equilibrium: a bracket by factors of four from
  ## b = 1, then a root-find in log b on its own cohort's R in its own field.
  R_of_b <- function(b) fitness(resident_at(b), resident_lma)
  R1 <- R_of_b(1)
  step <- if (R1 > 1) log(4) else -log(4)
  lb <- 0
  repeat {
    lb_next <- lb + step
    if (sign(log(R_of_b(exp(lb_next)))) != sign(log(R1))) break
    lb <- lb_next
    if (abs(lb) > log(4^8)) stop("no equilibrium birth rate within a factor 4^8 of 1")
  }
  w_root <- system.time(root <- uniroot(function(x) log(R_of_b(exp(x))), sort(c(lb, lb_next)),
                                        tol = 1e-10))[["elapsed"]]
  b_star <- exp(root$root)
  stamp("R(b = 1) = %.6g; b* = %.8g after %d evaluations in %.2f s, a cohort every %g yr",
        R1, b_star, root$iter, w_root, dt)

  w_res <- system.time(r <- resident_at(b_star))[["elapsed"]]
  sp <- r$scm$patch$species[[1]]
  stamp("resident at b*: run to %g yr in %.3f s; %d nodes, %d integrated, top-layer density %.4g",
        t_eq, w_res, sp$size, sp$size - sum(sp$retired), sp$top_layer$density)
  d <- do.call(rbind, lapply(lma_grid, function(lma) {
    w <- system.time(R <- fitness(r, lma))[["elapsed"]]
    stamp("mutant lma %.4f: R = %.6g in %.3f s", lma, R, w)
    data.frame(lma = lma, R = R, wall_s = w)
  }))
  d$R_resident <- d$R[d$lma == resident_lma]

  ## Inside one introduction interval, and a century either side.
  check_lma <- c(0.06, resident_lma, 0.1)
  fractions <- c(0, 0.25, 0.5, 0.75)
  ph <- sapply(fractions, function(f) {
    rr <- resident_at(b_star, t_eq + f * dt)
    sapply(check_lma, function(lma) fitness(rr, lma))
  })
  dimnames(ph) <- list(paste0("lma ", check_lma), paste0("+", fractions, " dt"))
  stamp("R in the field held at four points inside one introduction interval:")
  print(signif(ph, 8))
  stamp("relative spread inside the interval: %s",
        paste(signif(apply(ph, 1, function(x) diff(range(x)) / mean(x)), 3), collapse = ", "))
  drift <- sapply(c(-100, 0, 100), function(dt_read) {
    fitness(resident_at(b_star, t_eq + dt_read), resident_lma)
  })
  stamp("the resident's R read at %g, %g and %g yr: %s", t_eq - 100, t_eq, t_eq + 100,
        paste(signif(drift, 7), collapse = ", "))

  ## Against the same cohort born into the running stand.
  ref <- data.frame(lma = check_lma, t(sapply(check_lma, function(lma) in_stand(b_star, lma))))
  ref$R_field <- d$R[match(check_lma, d$lma)]
  ref$rel_diff <- ref$R_field / ref$R_mutant - 1
  stamp("against the full simultaneous run at b* (cohorts born at %g, run out %g yr):",
        t_eq + 0.5 * dt, life)
  print(format(ref, digits = 6), row.names = FALSE)
  write.csv(d, file.path(out_dir, paste0(prefix, "-equilibrium-mutants.csv")), row.names = FALSE)
  write.csv(ref, file.path(out_dir, paste0(prefix, "-equilibrium-check.csv")), row.names = FALSE)
}
stamp("done")
