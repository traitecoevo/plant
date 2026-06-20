## make_figure_gridartifact.R   (EXPERIMENTAL)
## ---------------------------------------------------------------------------
## "Forced to run" is not "valid": the PPA hard-step under fixed-step Euler.
##
## The adaptive solver REFUSES the hard-step shading variants (discontinuous
## light profile -> "Cannot achieve the desired accuracy"). A fixed-step
## forward-Euler integrator (control$fixed_time_step) will instead march
## straight across the discontinuity and return numbers -- but those numbers
## are an artefact of where the step grid falls, not a solution.
##
## This figure makes that concrete for PPA at two layer-smoothings:
##   (a) the SELECTION GRADIENT over resident LMA -- smooth (s=0.3) gives one
##       clean attractor; the hard step (s=0) is erratic and its whole sign
##       pattern CHANGES between a daily and a half-daily step.
##   (b) the underlying cause -- a single resident's seed output CONVERGES under
##       step refinement for the smooth profile but does NOT for the hard step.
##
## Selection gradients use the community / low-density invasion-fitness method
## (resident at demographic equilibrium + two mutants at tiny birth rate, one
## multi-species run_scm), which needs no RK45 mutant cache and so runs under
## fixed_time_step. Validated against scm$run_mutant to ~0.1% (adaptive).
##
##   Rscript vignettes/examples/DunFalster2026/make_figure_gridartifact.R
## ---------------------------------------------------------------------------

suppressMessages({
  devtools::load_all("/Users/z2209343/GitHub/packages/plant/plant-dev1")
  library(dplyr); library(ggplot2); library(patchwork)
})

here    <- "/Users/z2209343/GitHub/packages/plant/plant-dev1/vignettes/examples/DunFalster2026"
fig_dir <- file.path(here, "figures")       # git-ignored: outputs + data cache
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
out_ms  <- "/Users/z2209343/Library/CloudStorage/OneDrive-UNSW/atelier/2-research/active/Dun-competition/ms/curr/figures"

save_fig <- function(plot, basename, width, height) {
  dirs <- c(fig_dir, if (dir.exists(out_ms)) out_ms)
  for (d in dirs) {
    ggsave(file.path(d, paste0(basename, ".pdf")), plot, width = width,
           height = height, device = cairo_pdf)
    ggsave(file.path(d, paste0(basename, ".png")), plot, width = width,
           height = height, dpi = 300)
  }
  message("wrote ", basename, ".{pdf,png}")
}

theme_pub <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        legend.key.height = unit(0.9, "lines"))

p_fit <- scm_base_parameters("FF16"); p_fit$max_patch_lifetime <- 50

mk_ctrl <- function(dt, smoothing) {
  ctrl <- Control(); ctrl$shading_model <- "ppa"
  ctrl$ppa_layer_smoothing <- smoothing; ctrl$fixed_time_step <- dt
  ctrl
}
equilibrium_br <- function(ctrl, lma, n_iter = 6, br = 1) {
  for (i in seq_len(n_iter)) {
    res <- expand_parameters(trait_matrix(lma, "lma"), p_fit, birth_rate_list = br)
    br  <- run_scm(res, ctrl = ctrl)$offspring_production
  }
  br
}
## community / low-density selection gradient d log(R0_mutant)/d(LMA)
grad_community <- function(dt, smoothing, lma, delta = 0.0015, eps = 1e-4) {
  ctrl <- mk_ctrl(dt, smoothing)
  br   <- equilibrium_br(ctrl, lma)
  comm <- expand_parameters(
    trait_matrix(c(lma, lma - delta, lma + delta), "lma"), p_fit,
    birth_rate_list = c(br, eps, eps))
  op <- run_scm(comm, ctrl = ctrl)$offspring_production
  f  <- log(op[2:3] / eps)
  (f[2] - f[1]) / (2 * delta)
}

## ---- (a) gradient over LMA, three settings (cached) -----------------------
DAILY <- 1 / 365; HALF <- 1 / 730
settings <- tibble::tribble(
  ~setting,                ~dt,   ~smoothing,
  "PPA smoothed (daily)",  DAILY, 0.3,
  "PPA hard step (daily)", DAILY, 0,
  "PPA hard step (half-daily)", HALF, 0)
lma_grid <- seq(0.07, 0.12, length.out = 13)

cache_a <- file.path(fig_dir, "_gridartifact_gradient.csv")
if (file.exists(cache_a)) {
  grad_df <- read.csv(cache_a)
} else {
  grad_df <- bind_rows(lapply(seq_len(nrow(settings)), function(i) {
    s <- settings[i, ]; message("gradient: ", s$setting)
    g <- vapply(lma_grid, function(l)
      tryCatch(grad_community(s$dt, s$smoothing, l), error = function(e) NA_real_),
      numeric(1))
    data.frame(setting = s$setting, lma = lma_grid, gradient = g)
  }))
  write.csv(grad_df, cache_a, row.names = FALSE)
}
grad_df$setting <- factor(grad_df$setting, levels = settings$setting)

## ---- (b) single-resident seed output vs step size (cached) ----------------
cache_b <- file.path(fig_dir, "_gridartifact_convergence.csv")
if (file.exists(cache_b)) {
  conv_df <- read.csv(cache_b)
} else {
  dts <- 1 / c(46, 91, 182, 365, 730, 1460)
  op_at <- function(dt, smoothing) {
    res <- expand_parameters(trait_matrix(0.1, "lma"), p_fit, birth_rate_list = 1)
    tryCatch(run_scm(res, ctrl = mk_ctrl(dt, smoothing))$offspring_production,
             error = function(e) NA_real_)
  }
  conv_df <- bind_rows(lapply(c(0, 0.3), function(sm) {
    data.frame(smoothing = sm, steps_per_year = 1 / dts,
               offspring = vapply(dts, op_at, numeric(1), smoothing = sm))
  }))
  conv_df$variant <- ifelse(conv_df$smoothing == 0, "PPA hard step (s=0)",
                            "PPA smoothed (s=0.3)")
  write.csv(conv_df, cache_b, row.names = FALSE)
}

## ===========================================================================
cols <- c("PPA smoothed (daily)" = "#D55E00",
          "PPA hard step (daily)" = "#0072B2",
          "PPA hard step (half-daily)" = "#56B4E9")

p_a <- ggplot(grad_df, aes(lma, gradient, colour = setting)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  scale_colour_manual(values = cols, name = NULL) +
  coord_cartesian(ylim = c(-500, 500)) +
  labs(x = "Resident LMA", y = "Selection gradient  d log(R0)/d(LMA)",
       tag = "(a)",
       subtitle = "Smoothed PPA: one clean attractor. Hard step: erratic, and the\nsign pattern flips when the step is halved (clipped; |grad| up to ~9000).") +
  theme_pub + theme(legend.position = "bottom") +
  guides(colour = guide_legend(nrow = 1))

## free y-scales: the hard step's ~10x bounce (O(0.01-0.07)) is invisible next to
## the smoothed value (~25), so each variant needs its own scale to show the
## convergence (smoothed) vs non-convergence (hard) honestly.
p_b <- ggplot(conv_df, aes(steps_per_year, offspring, colour = variant)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  scale_x_log10() +
  facet_wrap(~variant, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = c("PPA hard step (s=0)" = "#0072B2",
                                 "PPA smoothed (s=0.3)" = "#D55E00"),
                      guide = "none") +
  labs(x = "Steps per year  (1 / fixed_time_step, log scale)",
       y = "Resident seed output", tag = "(b)",
       subtitle = "Smoothed converges under step refinement; the hard step does not (note each panel's own y-scale).") +
  theme_pub + theme(strip.text = element_text(face = "bold"))

fig <- (p_a / p_b) + plot_annotation(caption = paste(
  "Experimental. The PPA hard step (ppa_layer_smoothing = 0) is the literal field",
  "discretisation; the adaptive solver refuses it.\nForcing it to run with fixed-step",
  "forward Euler returns numbers that depend on the arbitrary step grid -- not a",
  "converged\nsolution. The smoothed variant (the runnable, reported version) is",
  "well-behaved. 'Does not run' understates it:\neven when forced to run, the hard",
  "step has no convergent evolutionary endpoint.")) &
  theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"))

save_fig(fig, "plant-grid-artifact", width = 8, height = 8)
message("\nDone.")
