## Mutants in the resident's equilibrium field (#519): a continuous patch with
## constant drivers, no disturbance and a cohort a year, run to its fixed point.
## There the resident's environment does not change, so a mutant needs no resident
## run of its own and no replay. Run the resident once to its fixed point, copy
## its environment (light field, and soil water for TF24), and integrate ONE
## mutant cohort in it until it is dead: its lifetime offspring is the mutant's
## fitness, because every arrival sees the same environment.
##
## Checked two ways: against the full simultaneous run, in which a cohort of the
## resident and of each checked mutant (at a birth rate 1e-6 of the resident's)
## is born at the reading year into the running stand and run out, so it sees the
## environment the stand actually makes; and the field read at four phases of the
## annual introduction cycle, which says how far the stand is from a strict
## fixed point between introductions.
##
## FF16 runs on the height coordinate by default: its equilibrium birth rate is
## near 70, and its stand overflows on the birth-date coordinate at b >= 10, a
## develop defect recorded in the plan. The top layer needs the birth-date
## coordinate, so on the height coordinate there is none; the canopy top edges up
## by centimetres a century there. TF24 runs on the birth-date coordinate (#590),
## with the top layer and constant rainfall: pulses never settle.
##
##   PLANT_519_MODEL    "FF16" (default) or "TF24"
##   PLANT_519_COORD    "height" or "birth_date" (FF16 only; TF24 is always birth_date)
##   PLANT_519_ABSORB   control$absorb_epsilon on the birth-date coordinate
##                      (0.02 for FF16, 0.05 for TF24)
##   PLANT_519_RAIN     TF24's constant rainfall, m per year (1)
##   PLANT_519_T_EQ     year at which the field is read (400)
##   PLANT_519_LIFE     years the reference run carries its checked cohorts (400)
##   PLANT_519_LMA      comma-separated mutant lma values
##   PLANT_519_RTOL     relative tolerance of the cohort integration (1e-8)
##   PLANT_519_OUT      output directory

suppressPackageStartupMessages(library(plant))
num_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(v)
}
vec_env <- function(name, default) {
  v <- Sys.getenv(name, "")
  if (v == "") default else as.numeric(strsplit(v, ",")[[1]])
}
model    <- Sys.getenv("PLANT_519_MODEL", "FF16")
tf24     <- identical(model, "TF24")
coord    <- if (tf24) "birth_date" else Sys.getenv("PLANT_519_COORD", "height")
birth_date <- identical(coord, "birth_date")
absorb   <- if (birth_date) num_env("PLANT_519_ABSORB", if (tf24) 0.05 else 0.02) else 0
rain     <- num_env("PLANT_519_RAIN", 1)
t_eq     <- num_env("PLANT_519_T_EQ", 400)
life     <- num_env("PLANT_519_LIFE", 400)
lma_grid <- vec_env("PLANT_519_LMA", c(0.04, 0.05, 0.06, 0.07, 0.0825, 0.1, 0.125, 0.15, 0.2))
rtol     <- num_env("PLANT_519_RTOL", 1e-8)
out_dir  <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}
stamp("plant %s; %s, resident lma 0.0825, %s coordinate, absorb_epsilon %g%s, field read at %g yr, rtol %g",
      as.character(packageVersion("plant")), model, coord, absorb,
      if (tf24) sprintf(", rainfall %g m/yr", rain) else "", t_eq, rtol)

resident_lma <- 0.0825
env_type <- if (tf24) "TF24_Env" else "FF16_Env"
## TF24's soil ODE is stiff for an explicit stepper without a cap on the step.
make_control <- function(absorb_epsilon) {
  if (tf24) {
    Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb_epsilon, ode_step_size_max = 0.05)
  } else {
    Control(node_density_in_birth_date = birth_date, absorb_epsilon = absorb_epsilon)
  }
}
ctrl <- make_control(absorb)
make_env <- function() {
  env <- Environment(model)
  if (tf24) env$extrinsic_drivers_set_constant("rainfall", rain)
  env
}
add_lma <- function(p, lma, b) {
  if (tf24) add_strategies(p, trait_matrix(lma, "lma"), hyperpar = TF24_hyperpar, birth_rate = b)
  else add_strategies(p, trait_matrix(lma, "lma"), birth_rate = b)
}

## Introductions every year from `phase`, so a run stopped at an integer time
## reads the stand `1 - phase` of the way through its annual cycle.
resident_parameters <- function(b, t_max, phase = 0) {
  p <- scm_base_parameters(model)
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- t_max
  p <- add_lma(p, resident_lma, b)
  p$node_schedule_times <- list(node_schedule_times_periodic(phase, t_max, 1))
  p
}

## ---- 1. the resident's environment at its fixed point -----------------------
resident_field <- function(b, phase = 0) {
  p <- resident_parameters(b, t_eq + 1, phase)
  scm <- plant:::SCM(model, env_type)(p, make_env(), plant:::empty_events(), ctrl)
  scm$run_to(t_eq)
  patch <- scm$patch
  patch$compute_environment()
  sp <- patch$species[[1]]
  list(env = patch$environment, time = patch$time, nodes = sp$size,
       live = sp$size - sum(sp$retired), top_density = sp$top_layer$density)
}

## ---- 2. one cohort in a frozen environment ----------------------------------
## The node's own ODE minus its log-density, which fitness does not read and
## whose growth-gradient term is a finite difference noisy enough to stall a
## stiff solver. odelia's Dormand--Prince over an R right-hand side, because a
## step that empties TF24's storage pool past zero has to be rejected and retried
## smaller, as plant's own stepper does; deSolve has no way to say so. Integrated
## in chunks until the cohort is dead to 1e-12.
strategy_of <- function(lma) {
  add_lma(resident_parameters(1, 1), lma, 1)$strategies[[2]]
}
cohort_fitness <- function(strategy, env) {
  node <- plant:::Node(model, env_type)(strategy)
  node$compute_initial_conditions(env, 1, 1)
  y_full <- node$ode_state
  keep <- which(node$ode_names != "log_density")
  names_kept <- node$ode_names[keep]
  i_storage <- match("storage", names_kept)
  calls <- 0L
  rhs <- function(t, y, parms) {
    if (!is.na(i_storage) && y[[i_storage]] < 0) odelia::domain_error("storage below zero")
    calls <<- calls + 1L
    y_full[keep] <- y
    node$ode_state <- y_full
    node$compute_rates(env, 1)
    node$ode_rates[keep]
  }
  y <- setNames(y_full[keep], names_kept)
  age <- 0
  repeat {
    out <- odelia::ode_solve(rhs, y, c(age, age + 100), rtol = rtol, atol = 1e-14)
    y <- setNames(out[2, -1], names_kept)
    age <- age + 100
    if (exp(-y[["mortality"]]) < 1e-12 || age >= 5000) break
  }
  c(R = y[["offspring_produced_survival_weighted"]] * strategy$pars$S_D,
    age = age, survival = exp(-y[["mortality"]]), height = y[["height"]], calls = calls)
}

## ---- 3. the reference: the same cohorts born into the running stand ---------
## The resident introduced every year as before, plus the mutant at the reading
## year at a birth rate 1e-6 of the resident's; the cohorts born then are run out
## for `life` years and scored by their own lifetime output. With the top layer
## on, a window over the reading year exempts them from absorption; without it
## nothing is absorbed.
simultaneous <- function(b, lma) {
  p <- resident_parameters(b, t_eq + life)
  p <- add_lma(p, lma, b * 1e-6)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, t_eq + life, 1), c(t_eq, t_eq + 1))
  if (absorb > 0) {
    p <- window_regime(p, t_eq, t_eq + 1)
    ctrl_ref <- ctrl
  } else {
    ctrl_ref <- make_control(0)
  }
  w <- system.time(scm <- run_scm(p, env = make_env(), ctrl = ctrl_ref))[["elapsed"]]
  sp <- scm$patch$species
  S_D <- p$strategies[[1]]$pars$S_D
  k <- which(sp[[1]]$node_times == t_eq)
  c(R_resident = sp[[1]]$net_reproduction_ratio_by_node[k] * S_D,
    R_mutant = sp[[2]]$net_reproduction_ratio_by_node[1] * p$strategies[[2]]$pars$S_D, wall = w)
}

## ---- 4. the resident's demographic equilibrium, by the field alone ----------
resident_strategy <- strategy_of(resident_lma)
R_of_b <- function(b) cohort_fitness(resident_strategy, resident_field(b)$env)[["R"]]
## A bracket by factors of four from b = 1, then a root-find in log b.
bracket <- c(0, 0)
R1 <- R_of_b(1)
step <- if (R1 > 1) log(4) else -log(4)
lb <- 0
repeat {
  lb_next <- lb + step
  if (sign(log(R_of_b(exp(lb_next)))) != sign(log(R1))) break
  lb <- lb_next
  if (abs(lb) > log(4^6)) stop("no equilibrium birth rate within a factor 4^6 of 1")
}
bracket <- sort(c(lb, lb_next))
w_root <- system.time(root <- uniroot(function(lb) log(R_of_b(exp(lb))), bracket,
                                      tol = 1e-10))[["elapsed"]]
b_star <- exp(root$root)
stamp("R(b = 1) = %.6g; b* = %.8g after %d evaluations in %.2f s, each a resident run to %g yr and one cohort",
      R1, b_star, root$iter, w_root, t_eq)

## ---- 5. the landscape at b*, every mutant in the one field ------------------
w_field <- system.time(f <- resident_field(b_star))[["elapsed"]]
stamp("field at b*: resident run to %g yr in %.3f s; %d nodes, %d integrated, top-layer density %.4g",
      f$time, w_field, f$nodes, f$live, f$top_density)
rows <- lapply(lma_grid, function(lma) {
  s <- strategy_of(lma)
  w <- system.time(r <- cohort_fitness(s, f$env))[["elapsed"]]
  stamp("mutant lma %.4f: R = %.6g in %.3f s (%d rate calls, run to age %g, height %.2f m)",
        lma, r[["R"]], w, as.integer(r[["calls"]]), r[["age"]], r[["height"]])
  data.frame(lma = lma, R = r[["R"]], wall_s = w, calls = r[["calls"]], age = r[["age"]],
             height_final = r[["height"]])
})
d <- do.call(rbind, rows)
d$R_resident <- d$R[d$lma == resident_lma]

## ---- 6. the phase of the annual cycle ---------------------------------------
check_lma <- c(0.06, resident_lma, 0.1)
phases <- c(0, 0.25, 0.5, 0.75)
ph <- sapply(phases, function(phi) {
  env <- resident_field(b_star, phi)$env
  sapply(check_lma, function(lma) cohort_fitness(strategy_of(lma), env)[["R"]])
})
dimnames(ph) <- list(paste0("lma ", check_lma), paste0("phase ", phases))
stamp("R read from the field at four phases of the introduction cycle:")
print(signif(ph, 8))
stamp("relative spread across phases: %s",
      paste(signif(apply(ph, 1, function(x) diff(range(x)) / mean(x)), 3), collapse = ", "))

## ---- 6b. drift: the same read a century either side ---------------------------
drift <- sapply(c(-100, 0, 100), function(dt) {
  t_read <- t_eq
  t_eq <<- t_read + dt
  on.exit(t_eq <<- t_read)
  cohort_fitness(resident_strategy, resident_field(b_star)$env)[["R"]]
})
stamp("the resident's R read at %g, %g and %g yr: %s", t_eq - 100, t_eq, t_eq + 100,
      paste(signif(drift, 7), collapse = ", "))

## ---- 7. against the full simultaneous run ------------------------------------
ref <- t(sapply(check_lma, function(lma) simultaneous(b_star, lma)))
ref <- data.frame(lma = check_lma, ref)
ref$R_field <- d$R[match(check_lma, d$lma)]
ref$rel_diff <- ref$R_field / ref$R_mutant - 1
ref$ratio_rel_diff <- (ref$R_field / d$R_resident[1]) / (ref$R_mutant / ref$R_resident) - 1
stamp("against the full simultaneous run at b* (cohorts born at %g, run out %g yr):", t_eq, life)
print(format(ref, digits = 6), row.names = FALSE)

prefix <- paste0(tolower(model), "-equilibrium-")
suffix <- if (!tf24 && birth_date) "-birth-date" else ""
write.csv(d, file.path(out_dir, paste0(prefix, "mutants", suffix, ".csv")), row.names = FALSE)
write.csv(ref, file.path(out_dir, paste0(prefix, "check", suffix, ".csv")), row.names = FALSE)
saveRDS(list(b_star = b_star, landscape = d, phases = ph, drift = drift, check = ref),
        file.path(out_dir, paste0(prefix, "mutants", suffix, ".rds")))
stamp("done")
