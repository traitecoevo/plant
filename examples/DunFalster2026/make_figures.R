## make_figures.R
## ---------------------------------------------------------------------------
## Publication figures for the "canopy approximations" manuscript
## (Dun & Falster 2026). All analysis is lifted verbatim from the plant
## vignette vignettes/canopy_methods.Rmd (chunks: lai, equil-table, fitness,
## compbox, softbox-fitness, sizedist, vprofile). The vignette is the source
## of truth; this script only makes the panels standalone and publication-grade.
##
## Run with the local dev build:
##   Rscript examples/DunFalster2026/make_figures.R
## ---------------------------------------------------------------------------

suppressMessages({
  devtools::load_all("/Users/z2209343/GitHub/packages/plant/plant-dev1")
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260620)

## ---------------------------------------------------------------------------
## Output locations
## ---------------------------------------------------------------------------
out_repo <- "/Users/z2209343/GitHub/packages/plant/plant-dev1/examples/DunFalster2026/figures"
out_ms   <- "/Users/z2209343/Library/CloudStorage/OneDrive-UNSW/atelier/2-research/active/Dun-competition/ms/curr/figures"
dir.create(out_repo, showWarnings = FALSE, recursive = TRUE)
dir.create(out_ms,   showWarnings = FALSE, recursive = TRUE)

## Save each figure as both PDF (vector; matches the filenames in the brief)
## and PNG (~300 dpi; matches the STYLE note). Written to repo + manuscript.
save_fig <- function(plot, basename, width, height) {
  for (dir in c(out_repo, out_ms)) {
    ggsave(file.path(dir, paste0(basename, ".pdf")), plot,
           width = width, height = height, device = cairo_pdf)
    ggsave(file.path(dir, paste0(basename, ".png")), plot,
           width = width, height = height, dpi = 300)
  }
  message("wrote ", basename, ".{pdf,png}")
}

## ---------------------------------------------------------------------------
## ONE place to set model ordering, colours and (swappable) labels
## ---------------------------------------------------------------------------
## Master ordering for every panel where shading model is an aesthetic.
## flat-top-soft-box sits between crown-centre and ppa so the ladder reads
## fine -> wrong (box) -> broken (ppa) left-to-right.
model_levels <- c("deep-crown", "mean-light", "crown-centre",
                  "flat-top-soft-box", "ppa", "flat-top-box")

## Colourblind-friendly (Okabe-Ito). deep-crown = black (the reference/default).
model_cols <- c(
  "deep-crown"        = "#000000",
  "mean-light"     = "#009E73",
  "crown-centre"          = "#0072B2",
  "ppa"               = "#D55E00",
  "flat-top-box"      = "#E69F00",
  "flat-top-soft-box" = "#CC79A7"
)

## Display labels: the manuscript naming scheme (paper <-> plant code). Changing
## this one vector relabels every axis, legend, facet and table below.
##
## Key point the names serve: plant's `crown-centre` keeps a SMOOTH (Yokozawa) shade
## and only flattens the *assimilation* calculation, so in the paper it is
## "Crown-Centre" -- NOT the field's "Flat Top", which means the box / thin-layer
## *shade* (plant's flat-top-box / flat-top-soft-box). Crown-Centre and Flat Top
## share the same crown-centre assimilation and differ only in the shade they
## cast, so the contrast isolates the shade cast, not the assimilation rule.
model_labels <- c(
  "deep-crown"        = "Deep Crown",
  "mean-light"     = "Mean-Light",
  "crown-centre"          = "Crown-Centre",
  "ppa"               = "PPA (smoothed)",
  "flat-top-box"      = "Flat Top (hard step)",
  "flat-top-soft-box" = "Flat Top (smoothed)"
)

## Two-colour scale for the fitness seedings (panel c).
seeding_cols <- c("same seed rain" = "#4477AA",
                  "each at its own equilibrium" = "#EE6677")

mk_factor <- function(x) factor(x, levels = model_levels)

scale_model_colour <- function(...) scale_colour_manual(
  values = model_cols, labels = model_labels, ..., drop = TRUE)
scale_model_fill <- function(...) scale_fill_manual(
  values = model_cols, labels = model_labels, ..., drop = TRUE)

theme_pub <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.key.height = unit(0.9, "lines"))

## ===========================================================================
## Shared setup: a single FF16 species, LMA = 0.0825 (exactly as the vignette)
## ===========================================================================
params <- scm_base_parameters("FF16")
patch  <- expand_parameters(trait_matrix(0.0825, "lma"), params)

## The four smooth-competition models (used for the supplementary panels and
## the bulk of the comparison) ...
smooth_models <- c("deep-crown", "mean-light", "crown-centre", "ppa")
## ... plus the soft-box, which runs but "goes wrong", added to the consequence
## ladder as the second-last column.
ladder_models <- c("deep-crown", "mean-light", "crown-centre",
                   "flat-top-soft-box", "ppa")
stand_models  <- union(smooth_models, ladder_models)

run_model <- function(model) {
  ctrl <- Control()
  ctrl$shading_model <- model
  run_scm(patch, ctrl = ctrl, collect = TRUE)
}

message("running the stands ...")
results <- lapply(stand_models, run_model)
names(results) <- stand_models

species_all <- dplyr::bind_rows(
  lapply(stand_models, function(m) dplyr::mutate(results[[m]]$species, method = m))
)
species_all$method <- mk_factor(species_all$method)

## Light-extinction coefficient: E = exp(-k_I * LAI), so LAI = -log(E)/k_I.
k_I <- results[[1]]$p$strategies[[1]]$k_I

## ===========================================================================
## FIGURE 1 — plant-consequence-ladder  (fine -> diverging -> broken)
## ===========================================================================

## ---- (a) Stand dynamics: leaf area index through patch age  [chunk: lai] ----
lai_all <- dplyr::bind_rows(
  lapply(ladder_models, function(m) {
    results[[m]]$env$light_availability |>
      dplyr::group_by(time) |>
      dplyr::slice_min(height, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::transmute(method = m, time, lai = -log(light_availability) / k_I)
  })
)
lai_all$method <- mk_factor(lai_all$method)

p_lai <- ggplot(lai_all, aes(time, lai, colour = method)) +
  geom_line(linewidth = 0.9) +
  scale_model_colour(name = "Shading model") +
  labs(x = "Patch age (years)", y = "Leaf area index",
       tag = "(a)", subtitle = "Stand dynamics: fine") +
  theme_pub

## ---- (b) Demographic equilibrium: equilibrium seed rain  [chunk: equil-table]
## Hardcoded equilibria from the vignette (slow models are expensive to re-derive).
equilibrium_birth_rate <- c(
  "deep-crown"        = 1.662,
  "mean-light"     = 2.178,
  "crown-centre"          = 2.185,
  "ppa"               = 18.835,
  "flat-top-soft-box" = 3.101
)

equil_df <- data.frame(
  method = factor(ladder_models, levels = model_levels),
  seed_rain = equilibrium_birth_rate[ladder_models]
)

## Lollipop on a log scale makes the order-of-magnitude spread obvious.
p_equil <- ggplot(equil_df, aes(method, seed_rain, colour = method)) +
  geom_segment(aes(xend = method, y = min(seed_rain) * 0.6, yend = seed_rain),
               linewidth = 0.9) +
  geom_point(size = 3.5) +
  geom_text(aes(label = formatC(seed_rain, format = "f", digits = 2)),
            vjust = -1.1, size = 3, show.legend = FALSE) +
  scale_model_colour(name = "Shading model") +
  scale_y_log10(expand = expansion(mult = c(0.02, 0.18))) +
  scale_x_discrete(labels = model_labels) +
  labs(x = NULL, y = "Equilibrium seed rain\n(log scale)",
       tag = "(b)", subtitle = "Demographic equilibrium: diverging") +
  theme_pub +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

## ---- (c) Invasion fitness landscapes, two seedings  [chunk: fitness] --------
p_fit <- scm_base_parameters("FF16")
p_fit$max_patch_lifetime <- 50
lma_resident <- 0.0825
mutant_lma <- seq(0.07, 0.11, length.out = 60)

fitness_landscape <- function(model, birth_rate) {
  ctrl <- Control()
  ctrl$shading_model <- model
  ctrl$save_RK45_cache <- TRUE                      # required by run_mutant()
  resident <- expand_parameters(trait_matrix(lma_resident, "lma"), p_fit,
                                birth_rate_list = birth_rate)
  scm <- run_scm(resident, ctrl = ctrl)             # caches the light environment
  mutants <- expand_parameters(trait_matrix(mutant_lma, "lma"), resident,
                               birth_rate_list = rep(1, length(mutant_lma)))
  scm$run_mutant(mutants)
  r0 <- scm$net_reproduction_ratios                 # resident first, then mutants
  data.frame(method = model, lma = mutant_lma,
             fitness = log(tail(r0, length(mutant_lma))))
}

message("computing fitness landscapes (same seed rain) ...")
common_rain <- equilibrium_birth_rate[["deep-crown"]]
fit_same  <- dplyr::bind_rows(lapply(ladder_models, fitness_landscape,
                                     birth_rate = common_rain))
message("computing fitness landscapes (each at its own equilibrium) ...")
fit_equil <- dplyr::bind_rows(lapply(ladder_models, function(m)
  fitness_landscape(m, equilibrium_birth_rate[[m]])))

fit_two <- dplyr::bind_rows(
  dplyr::mutate(fit_same,  seeding = "same seed rain"),
  dplyr::mutate(fit_equil, seeding = "each at its own equilibrium")
)
fit_two$method  <- mk_factor(fit_two$method)
fit_two$seeding <- factor(fit_two$seeding,
                          levels = c("same seed rain", "each at its own equilibrium"))

## The columns shown are the runnable models. Note in the relevant facets that
## their hard/literal siblings (the box step, and an unsmoothed PPA) fail to run.
ladder_fail <- data.frame(
  method = factor(c("flat-top-soft-box", "ppa"), levels = model_levels),
  lma = 0.090,
  fitness = -13,
  label = c("Flat Top (hard step)\nfails to run",
            "PPA (hard step)\nfails to run")
)

p_fitness <- ggplot(fit_two, aes(lma, fitness, colour = seeding)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
  geom_line(linewidth = 0.8) +
  geom_text(data = ladder_fail, inherit.aes = FALSE,
            aes(lma, fitness, label = label),
            size = 2.5, colour = "grey35", lineheight = 0.9) +
  facet_wrap(~method, nrow = 1, labeller = as_labeller(model_labels)) +
  scale_colour_manual(values = seeding_cols, name = NULL) +
  labs(x = "Mutant LMA", y = "Invasion fitness, log(R0)",
       tag = "(c)", subtitle = "Invasion fitness: broken (Flat Top smoothed biased, PPA jagged)") +
  theme_pub +
  theme(legend.position = "bottom")

ladder <- (p_lai / p_equil / p_fitness) +
  plot_layout(heights = c(1, 1, 1.05)) +
  plot_annotation(caption = paste(
    "Columns are the runnable models. The literal field-standard discretisations —",
    "Flat Top (hard step) and PPA (hard step) — fail to run\n(discontinuous light profile)",
    "and cannot be shown; the Flat Top (smoothed) and PPA columns use the default layer smoothing.")) &
  theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"))

save_fig(ladder, "plant-consequence-ladder", width = 8.6, height = 10.5)

## ===========================================================================
## FIGURE 2 — competition-profile-and-fitness  (the mechanism, side by side)
## ===========================================================================

## ---- (left) shade cast vs height  [chunk: compbox] --------------------------
comp_profile <- function(model, height = 10) {
  s <- FF16_Strategy()
  s$control$shading_model <- model
  ind <- FF16_Individual(s)
  ind$set_state("height", height)
  z <- seq(0, height, length.out = 400)
  data.frame(model = model, z = z,
             competition = sapply(z, function(zz) ind$compute_competition(zz)))
}
comp <- dplyr::bind_rows(
  comp_profile("crown-centre"),
  comp_profile("flat-top-box"),
  comp_profile("flat-top-soft-box")
)
comp$model <- mk_factor(comp$model)

p_comp <- ggplot(comp, aes(competition, z, colour = model)) +
  geom_path(linewidth = 0.9) +
  scale_model_colour(name = "Model") +
  labs(x = "Shade cast (competition) by a 10 m plant", y = "Height (m)",
       tag = "(a)", subtitle = "Competition profile") +
  theme_pub +
  theme(legend.position = "bottom", legend.direction = "vertical")

## ---- (right) fitness landscapes, each at its own equilibrium [softbox-fitness]
box_models <- c("deep-crown", "crown-centre", "flat-top-soft-box")
message("computing soft-box fitness comparison ...")
fit_compare <- dplyr::bind_rows(lapply(box_models, function(m)
  fitness_landscape(m, equilibrium_birth_rate[[m]])))
fit_compare$method <- mk_factor(fit_compare$method)

## Flat Top (hard step) casts the box shade shown in panel (a) but yields no
## fitness curve here: summed over cohorts its step makes the patch light
## profile discontinuous, so no stand can be solved. Flag that absence.
softbox_fail <- data.frame(
  lma = min(fit_compare$lma),
  fitness = max(fit_compare$fitness, na.rm = TRUE),
  label = "✗  Flat Top (hard step): no curve —\nstep shade (panel a) gives a discontinuous\nlight profile, so no stand can be solved")

p_softbox <- ggplot(fit_compare, aes(lma, fitness, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 1) +
  geom_text(data = softbox_fail, inherit.aes = FALSE,
            aes(lma, fitness, label = label), hjust = 0, vjust = 1,
            size = 2.6, colour = model_cols[["flat-top-box"]], lineheight = 0.95) +
  scale_model_colour(name = "Model") +
  labs(x = "Mutant LMA", y = "Invasion fitness, log(R0)",
       tag = "(b)", subtitle = "Invasion fitness at each own equilibrium") +
  theme_pub +
  theme(legend.position = "bottom", legend.direction = "vertical")

mechanism <- (p_comp | p_softbox) + plot_layout(widths = c(1, 1.1))
save_fig(mechanism, "competition-profile-and-fitness", width = 11, height = 5)

## ===========================================================================
## SUPPLEMENTARY FIGURE 3 — plant-size-distribution  [chunk: sizedist]
## ===========================================================================
species_smooth <- species_all |>
  dplyr::filter(method %in% smooth_models) |>
  dplyr::mutate(method = droplevels(method))
p_size <- plot_size_distribution(species_smooth) +
  facet_wrap(~method, nrow = 1, labeller = as_labeller(model_labels)) +
  theme_pub +
  theme(legend.position = "bottom")

save_fig(p_size, "plant-size-distribution", width = 10, height = 3.6)

## ===========================================================================
## SUPPLEMENTARY FIGURE 4 — plant-vertical-light-profile  [chunk: vprofile]
## ===========================================================================
d <- Control()$ppa_layer_optical_depth   # layer thickness (optical depth)
w <- Control()$ppa_layer_smoothing       # boundary smoothing fraction

## C1-smoothed staircase, mirroring FF16_Environment::smooth_floor()/step_light().
smooth_floor <- function(u, w) {
  n <- floor(u)
  f <- u - n
  t <- pmax(0, f - (1 - w)) / w
  n + ifelse(f <= 1 - w, 0, t^2 * (3 - 2 * t))
}
ppa_step_light <- function(E, d, w) {
  ifelse(E >= 1 | E <= 0, E,
         exp(-d * smooth_floor(-log(E) / d, w)))
}

target_ages <- c(5, 10, 20, 60)

profile_all <- dplyr::bind_rows(
  lapply(smooth_models, function(m) {
    light <- results[[m]]$env$light_availability
    if (m == "ppa") {
      light$light_availability <- ppa_step_light(light$light_availability, d, w)
    }
    steps <- results[[m]]$steps
    keep <- sapply(target_ages, function(a) steps$step[which.min(abs(steps$time - a))])
    light |>
      dplyr::filter(step %in% keep) |>
      dplyr::mutate(method = m,
                    age = target_ages[match(step, keep)],
                    lai_above = -log(light_availability) / k_I)
  })
)
profile_all$method <- mk_factor(profile_all$method)

p_profile <- ggplot(profile_all,
                    aes(lai_above, height, colour = factor(age), group = age)) +
  geom_path(linewidth = 0.8) +
  facet_wrap(~method, nrow = 1, labeller = as_labeller(model_labels)) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = expression(LAI[above]~(m^2/m^2)), y = "Height (m)",
       colour = "Patch age\n(years)") +
  theme_pub

save_fig(p_profile, "plant-vertical-light-profile", width = 10, height = 3.6)

## ===========================================================================
## SUPPLEMENTARY — code <-> paper model-name mapping (for reproducibility)
## ===========================================================================
mapping <- data.frame(
  `Paper name` = c("Deep Crown", "PPA (smoothed)", "PPA (hard step)",
                   "Flat Top", "Mean-Light", "Crown-Centre"),
  `Meaning` = c(
    "smooth (Yokozawa) shade; assimilation integrated over crown depth",
    "layered, stepped light profile, boundaries smoothed (runnable; all PPA results)",
    "literal field discretisation; discontinuous, does not run",
    "box / thin-layer shade (field's usual sense; = Stage-1 MAESPA \"Flat Top\")",
    "smooth shade; light averaged over crown depth, one assimilation evaluation",
    "smooth shade; one assimilation evaluation at the crown centre"),
  `plant shading_model code` = c(
    "deep-crown",
    "ppa, ppa_layer_smoothing > 0 (default 0.3)",
    "ppa, ppa_layer_smoothing = 0",
    "flat-top-box (hard step), flat-top-soft-box (smoothed)",
    "mean-light  (plant-only, new)", "crown-centre  (plant-only, new)"),
  check.names = FALSE, stringsAsFactors = FALSE
)

## Reproducible markdown table (lives with the script in the repo).
md <- c(
  "# Canopy-model naming: paper <-> plant code", "",
  "| Paper name | Meaning | `plant` `shading_model` code |",
  "|---|---|---|",
  apply(mapping, 1, function(r)
    sprintf("| %s | %s | %s |", r[[1]], r[[2]], r[[3]])),
  "",
  "- The code label `crown-centre` matches the paper name **Crown-Centre** (smooth shade, crown-centre assimilation) — it is **not** the field's \"Flat Top\", which is the box / thin-layer shade (plant's `flat-top-box` / `flat-top-soft-box`).",
  "- **Flat Top** = the box family: \"Flat Top (hard step)\" (does not run) vs \"Flat Top (smoothed)\" (runs but biased).",
  "- **PPA** is one shading model (`ppa`) at two settings of `ppa_layer_smoothing`, not two `shading_model` strings: \"PPA (smoothed)\" (`> 0`, default 0.3 — runnable; all reported PPA results) vs \"PPA (hard step)\" (`= 0` — the literal field discretisation; discontinuous, does not run in the adaptive solver).",
  "- **Mean-Light** and **Crown-Centre** appear only in the dynamic (Stage 2) results.",
  "- Crown-Centre and Flat Top use the *same* crown-centre assimilation and differ only in the shade cast — so Crown-Centre giving usable fitness while Flat Top fails isolates the shade cast, not the assimilation rule."
)
writeLines(md, file.path(dirname(out_repo), "model-name-mapping.md"))

## Rendered table figure for the supplement.
tbl_theme <- gridExtra::ttheme_minimal(
  base_size = 9, padding = unit(c(5, 5), "mm"),
  core = list(fg_params = list(hjust = 0, x = 0.02)),
  colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold")))
tbl_grob <- gridExtra::tableGrob(mapping, rows = NULL, theme = tbl_theme)
for (dir in c(out_repo, out_ms)) {
  ggsave(file.path(dir, "plant-model-name-mapping.pdf"), tbl_grob,
         width = 10, height = 2.6, device = cairo_pdf)
  ggsave(file.path(dir, "plant-model-name-mapping.png"), tbl_grob,
         width = 10, height = 2.6, dpi = 300)
}
message("wrote plant-model-name-mapping.{pdf,png} and model-name-mapping.md")

message("\nAll figures written to:\n  ", out_repo, "\n  ", out_ms)
