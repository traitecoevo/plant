## Mutants in the resident's equilibrium field (#519): at a fixed point the
## resident's environment does not change, so a mutant needs no resident run of
## its own and no replay. Run the resident once to its fixed point, copy the
## light field, and integrate ONE mutant cohort in that field until it is dead:
## its lifetime offspring is the mutant's fitness, because every arrival sees the
## same field.
##
## Checked two ways: against the full simultaneous run, in which a cohort of the
## resident and of each checked mutant (at a birth rate 1e-6 of the resident's)
## is born at the reading year into the running stand and run out, so it sees the
## field the stand actually makes; and the field read at four phases of the
## annual introduction cycle, which says how far the stand is from a strict
## fixed point between introductions.
##
## The height coordinate by default: FF16's equilibrium birth rate is near 70,
## and this stand overflows on the birth-date coordinate at b >= 10, a develop
## defect recorded in the plan. The top layer needs the birth-date coordinate,
## so on the height coordinate there is none; the canopy top edges up by
## centimetres a century there, at densities that shade nothing.
##
##   PLANT_519_COORD    "height" (default) or "birth_date"
##   PLANT_519_ABSORB   control$absorb_epsilon on the birth-date coordinate (0.02)
##   PLANT_519_T_EQ     year at which the field is read (400)
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
coord    <- Sys.getenv("PLANT_519_COORD", "height")
birth_date <- identical(coord, "birth_date")
absorb   <- if (birth_date) num_env("PLANT_519_ABSORB", 0.02) else 0
t_eq     <- num_env("PLANT_519_T_EQ", 400)
lma_grid <- vec_env("PLANT_519_LMA", c(0.04, 0.05, 0.06, 0.07, 0.0825, 0.1, 0.125, 0.15, 0.2))
rtol     <- num_env("PLANT_519_RTOL", 1e-8)
out_dir  <- Sys.getenv("PLANT_519_OUT", "notes/scripts/continuous-patch-519/output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n", sep = " ")
  flush(stdout())
}
stamp("plant %s; FF16, resident lma 0.0825, %s coordinate, absorb_epsilon %g, field read at %g yr, rtol %g",
      as.character(packageVersion("plant")), coord, absorb, t_eq, rtol)

resident_lma <- 0.0825
ctrl <- Control(node_density_in_birth_date = birth_date, absorb_epsilon = absorb)
## The reference run never absorbs, so the checked cohorts keep their whole output.
ctrl_ref <- Control(node_density_in_birth_date = birth_date)

## Introductions every year from `phase`, so a run stopped at an integer time
## reads the stand `1 - phase` of the way through its annual cycle.
resident_parameters <- function(b, t_max, phase = 0) {
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- t_max
  p <- add_strategies(p, trait_matrix(resident_lma, "lma"), birth_rate = b)
  p$node_schedule_times <- list(node_schedule_times_periodic(phase, t_max, 1))
  p
}

## ---- 1. the resident's field at its fixed point -----------------------------
resident_field <- function(b, phase = 0) {
  p <- resident_parameters(b, t_eq + 1, phase)
  scm <- plant:::SCM("FF16", "FF16_Env")(p, Environment("FF16"), plant:::empty_events(), ctrl)
  scm$run_to(t_eq)
  patch <- scm$patch
  patch$compute_environment()
  sp <- patch$species[[1]]
  list(env = patch$environment, time = patch$time, nodes = sp$size,
       live = sp$size - sum(sp$retired), top_density = sp$top_layer$density)
}

## ---- 2. one cohort in a frozen field ----------------------------------------
## The node's own ODE minus its log-density, which fitness does not read and
## whose growth-gradient term is a finite difference noisy enough to stall a
## stiff solver. Integrated in chunks until the cohort is dead to 1e-12.
strategy_of <- function(lma) {
  p <- resident_parameters(1, 1)
  add_strategies(p, trait_matrix(lma, "lma"), birth_rate = 1)$strategies[[2]]
}
cohort_fitness <- function(strategy, env) {
  node <- plant:::Node("FF16", "FF16_Env")(strategy)
  node$compute_initial_conditions(env, 1, 1)
  y_full <- node$ode_state
  keep <- which(node$ode_names != "log_density")
  calls <- 0L
  rhs <- function(t, y, parms) {
    calls <<- calls + 1L
    y_full[keep] <- y
    node$ode_state <- y_full
    node$compute_rates(env, 1)
    list(node$ode_rates[keep])
  }
  y <- setNames(y_full[keep], node$ode_names[keep])
  age <- 0
  repeat {
    out <- deSolve::ode(y, c(age, age + 100), rhs, parms = NULL, method = "ode45",
                        rtol = rtol, atol = 1e-14, maxsteps = 1e6)
    y <- out[2, -1]
    age <- age + 100
    if (exp(-y[["mortality"]]) < 1e-12 || age >= 5000) break
  }
  c(R = y[["offspring_produced_survival_weighted"]] * strategy$pars$S_D,
    age = age, survival = exp(-y[["mortality"]]), height = y[["height"]], calls = calls)
}

## ---- 3. the reference: the same cohorts born into the running stand ---------
## The resident introduced every year as before, plus one cohort of the mutant at
## the reading year at a birth rate 1e-6 of the resident's; both cohorts born
## then are run out for `life` years and scored by their own lifetime output.
simultaneous <- function(b, lma, life = 400) {
  p <- resident_parameters(b, t_eq + life)
  p <- add_strategies(p, trait_matrix(lma, "lma"), birth_rate = b * 1e-6)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, t_eq + life, 1), t_eq)
  w <- system.time(scm <- run_scm(p, ctrl = ctrl_ref))[["elapsed"]]
  sp <- scm$patch$species
  S_D <- p$strategies[[1]]$pars$S_D
  k <- which(sp[[1]]$node_times == t_eq)
  c(R_resident = sp[[1]]$net_reproduction_ratio_by_node[k] * S_D,
    R_mutant = sp[[2]]$net_reproduction_ratio_by_node[1] * p$strategies[[2]]$pars$S_D, wall = w)
}

## ---- 4. the resident's demographic equilibrium, by the field alone ----------
resident_strategy <- strategy_of(resident_lma)
R_of_b <- function(b) cohort_fitness(resident_strategy, resident_field(b)$env)[["R"]]
w_root <- system.time(root <- uniroot(function(lb) log(R_of_b(exp(lb))), c(log(10), log(500)),
                                      tol = 1e-10))[["elapsed"]]
b_star <- exp(root$root)
stamp("b* = %.8g after %d evaluations in %.2f s, each a resident run to %g yr and one cohort",
      b_star, root$iter, w_root, t_eq)

## ---- 5. the landscape at b*, every mutant in the one field ------------------
w_field <- system.time(f <- resident_field(b_star))[["elapsed"]]
stamp("field at b*: resident run to %g yr in %.3f s; %d nodes, %d integrated, top-layer density %.4g",
      f$time, w_field, f$nodes, f$live, f$top_density)
rows <- lapply(lma_grid, function(lma) {
  s <- strategy_of(lma)
  w <- system.time(r <- cohort_fitness(s, f$env))[["elapsed"]]
  stamp("mutant lma %.4f: R = %.6f in %.3f s (%d rate calls, run to age %g, height %.2f m)",
        lma, r[["R"]], w, as.integer(r[["calls"]]), r[["age"]], r[["height"]])
  data.frame(lma = lma, R = r[["R"]], wall_s = w, calls = r[["calls"]], age = r[["age"]],
             height_final = r[["height"]])
})
d <- do.call(rbind, rows)
d$R_resident <- d$R[d$lma == resident_lma]

## ---- 6. the phase of the annual cycle ---------------------------------------
phases <- c(0, 0.25, 0.5, 0.75)
ph <- sapply(phases, function(phi) {
  env <- resident_field(b_star, phi)$env
  sapply(c(0.06, resident_lma, 0.1), function(lma) cohort_fitness(strategy_of(lma), env)[["R"]])
})
dimnames(ph) <- list(paste0("lma ", c(0.06, resident_lma, 0.1)), paste0("phase ", phases))
stamp("R read from the field at four phases of the introduction cycle:")
print(signif(ph, 8))
stamp("relative spread across phases: %s",
      paste(signif(apply(ph, 1, function(x) diff(range(x)) / mean(x)), 3), collapse = ", "))

## ---- 7. against the full simultaneous run ------------------------------------
check_lma <- c(0.06, resident_lma, 0.1)
ref <- t(sapply(check_lma, function(lma) simultaneous(b_star, lma)))
ref <- data.frame(lma = check_lma, ref)
ref$R_field <- d$R[match(check_lma, d$lma)]
ref$rel_diff <- ref$R_field / ref$R_mutant - 1
ref$ratio_rel_diff <- (ref$R_field / d$R_resident[1]) / (ref$R_mutant / ref$R_resident) - 1
stamp("against the full simultaneous run at b* (cohorts born at %g, run out 400 yr):", t_eq)
print(format(ref, digits = 6), row.names = FALSE)

suffix <- if (birth_date) "-birth-date" else ""
write.csv(d, file.path(out_dir, paste0("ff16-equilibrium-mutants", suffix, ".csv")), row.names = FALSE)
write.csv(ref, file.path(out_dir, paste0("ff16-equilibrium-check", suffix, ".csv")), row.names = FALSE)
saveRDS(list(b_star = b_star, landscape = d, phases = ph, check = ref),
        file.path(out_dir, paste0("ff16-equilibrium-mutants", suffix, ".rds")))
stamp("done")
