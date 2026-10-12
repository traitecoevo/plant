## The absorbing top layer on a continuous FF16 stand (#519): exact conservation
## at an absorption instant, and the stand under the layer at three heights
## against the unbounded stand. The table in the plan's section 7; its survival
## floor row came from code since removed (it is in this branch's history).
suppressPackageStartupMessages(library(plant))
trap <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
mk <- function(T, b = 1) { p <- scm_base_parameters("FF16"); p$patch_type <- "fixed"; p$max_patch_lifetime <- T
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = b); p$node_schedule_times <- list(node_schedule_times_periodic(0, T, 1)); p }
ctrl <- function(absorb = 0) Control(node_density_in_birth_date = TRUE, absorb_epsilon = absorb)

## 1. Exact conservation at the instant: a stand grown without absorption, so
## many cohorts stand above the absorbing height, loaded with it on.
burn <- plant:::SCM("FF16", "FF16_Env")(mk(400), Environment("FF16"), plant:::empty_events(), ctrl())
burn$run_to(300); state <- export_patch_state(burn)
state$ode_state <- c(state$ode_state, rep(0, 5 + 5))   # an empty top block (FF16: 5 states + 5)
p2 <- set_initial_state(mk(400), state)
patch <- plant:::Patch("FF16", "FF16_Env")(p2, Environment("FF16"), ctrl(absorb = 0.05))
total <- function(pt) { sp <- pt$species[[1]]; t <- c(sp$node_times, pt$time)
  trap(t, c(exp(sp$log_densities_state), exp(sp$new_node$log_density))) + sp$top_layer$density }
h_abs <- FF16_Strategy()$pars$hmat * (1 + log(19) / FF16_Strategy()$pars$a_f2)
sp <- patch$species[[1]]
cat(sprintf("before: %d cohorts, %d at or above h_abs = %.4f m, top density %g, total density %.15g\n", sp$size, sum(sp$heights >= h_abs), h_abs, sp$top_layer$density, total(patch)))
D0 <- total(patch)
patch$introduce_new_node(1L, patch$time)
sp <- patch$species[[1]]
cat(sprintf("after:  %d cohorts, %d retired, top density %.6g at height %.6f, total density %.15g, rel change %.2e\n", sp$size, sum(sp$retired), sp$top_layer$density, sp$top_layer$height, total(patch), abs(total(patch) / D0 - 1)))

## 2. The continuous stand, three boundaries.
cfgs <- list(list(lab = "no boundary", a = 0))
for (eps in c(0.05, 0.02, 0.01)) cfgs[[length(cfgs) + 1]] <- list(lab = sprintf("top layer, epsilon %g", eps), a = eps)
for (cfg in cfgs) {
  s <- plant:::SCM("FF16", "FF16_Env")(mk(1000), Environment("FF16"), plant:::empty_events(), ctrl(cfg$a))
  s$snapshot_interval <- 50
  w <- system.time(s$run())[["elapsed"]]
  sp <- s$patch$species[[1]]; sn <- tidy_snapshots(s)$species
  late <- sn[sn$time >= 500, ]
  cat(sprintf("%-26s %.2f s; integrated %4d; top density %.4g at %s; from 500 yr: tallest %.5f to %.5f m, density %.7g to %.7g; leaf area %.9g, biomass %.9g, R %.9g\n",
    cfg$lab, w, sp$size - sum(sp$retired), sp$top_layer$density, format(sp$top_layer$height, digits = 6),
    min(late$height_max), max(late$height_max), min(late$density), max(late$density),
    stand_census(s)[["leaf_area"]], stand_census(s)[["mass_above_ground"]], s$net_reproduction_ratios))
}
