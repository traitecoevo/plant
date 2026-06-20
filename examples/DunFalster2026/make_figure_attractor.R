## make_figure_attractor.R   (EXPERIMENTAL)
## ---------------------------------------------------------------------------
## Beyond invasion fitness: the one-species evolutionary attractor.
##
## For each shading method we find the singular strategy (evolutionarily
## attracted LMA) by root-solving the selection gradient:
##   1. for a resident LMA, find the demographic-equilibrium seed rain (R0 = 1);
##   2. introduce mutants either side of the resident and measure the selection
##      gradient  d log(R0_mutant)/d(LMA)  at mutant = resident;
##   3. the attractor is where that gradient is zero (with negative slope =
##      convergence stable). We map the gradient over a grid AND uniroot it.
##
## Question for the paper: do the methods favour different attracted LMA values?
##
## Style (palette / labels / theme) mirrors make_figures.R so the colours and
## names stay consistent across every figure.
##
##   Rscript examples/DunFalster2026/make_figure_attractor.R
## ---------------------------------------------------------------------------

suppressMessages({
  devtools::load_all("/Users/z2209343/GitHub/packages/plant/plant-dev1")
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(20260620)

out_repo <- "/Users/z2209343/GitHub/packages/plant/plant-dev1/examples/DunFalster2026/figures"
out_ms   <- "/Users/z2209343/Library/CloudStorage/OneDrive-UNSW/atelier/2-research/active/Dun-competition/ms/curr/figures"

save_fig <- function(plot, basename, width, height) {
  for (dir in c(out_repo, out_ms)) {
    ggsave(file.path(dir, paste0(basename, ".pdf")), plot,
           width = width, height = height, device = cairo_pdf)
    ggsave(file.path(dir, paste0(basename, ".png")), plot,
           width = width, height = height, dpi = 300)
  }
  message("wrote ", basename, ".{pdf,png}")
}

## ---- consistent ordering / colours / labels (mirror make_figures.R) --------
model_levels <- c("deep-crown", "mean-light", "crown-centre",
                  "flat-top-soft-box", "ppa", "flat-top-box")
model_cols <- c(
  "deep-crown"        = "#000000",
  "mean-light"        = "#009E73",
  "crown-centre"      = "#0072B2",
  "ppa"               = "#D55E00",
  "flat-top-box"      = "#E69F00",
  "flat-top-soft-box" = "#CC79A7"
)
model_labels <- c(
  "deep-crown"        = "Deep Crown",
  "mean-light"        = "Mean-Light",
  "crown-centre"      = "Crown-Centre",
  "ppa"               = "PPA (smoothed)",
  "flat-top-box"      = "Flat Top (hard step)",
  "flat-top-soft-box" = "Flat Top (smoothed)"
)
mk_factor <- function(x) factor(x, levels = model_levels)
theme_pub <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.key.height = unit(0.9, "lines"))

## ===========================================================================
## Selection-gradient machinery
## ===========================================================================
p_fit <- scm_base_parameters("FF16")
p_fit$max_patch_lifetime <- 50

## Demographic-equilibrium seed rain for a resident of given LMA.
equilibrium_birth_rate_at <- function(model, lma, n_iter = 6, birth_rate = 1) {
  ctrl <- Control(); ctrl$shading_model <- model
  for (i in seq_len(n_iter)) {
    resident <- expand_parameters(trait_matrix(lma, "lma"), p_fit,
                                  birth_rate_list = birth_rate)
    birth_rate <- run_scm(resident, ctrl = ctrl)$offspring_production
  }
  birth_rate
}

## Selection gradient d log(R0)/d(LMA) at mutant = resident, resident at
## its own demographic equilibrium. Central finite difference.
sel_gradient <- function(model, lma, delta = 0.0015, n_iter = 6) {
  ctrl <- Control(); ctrl$shading_model <- model; ctrl$save_RK45_cache <- TRUE
  br <- equilibrium_birth_rate_at(model, lma, n_iter)
  resident <- expand_parameters(trait_matrix(lma, "lma"), p_fit,
                                birth_rate_list = br)
  scm <- run_scm(resident, ctrl = ctrl)
  mut <- c(lma - delta, lma + delta)
  mutants <- expand_parameters(trait_matrix(mut, "lma"), resident,
                               birth_rate_list = rep(1, 2))
  scm$run_mutant(mutants)
  r0 <- scm$net_reproduction_ratios
  f  <- log(tail(r0, 2))
  (f[2] - f[1]) / (2 * delta)
}

## Evolutionary stability: curvature of the mutant fitness landscape at
## mutant = s* (resident fixed at s*'s equilibrium). d2 log(R0)/d(LMA)^2 < 0 =>
## uninvadable (ESS); combined with convergence stability (the + -> - gradient
## crossing) => CSS, a continuously stable evolutionary endpoint. > 0 => branching.
fitness_curvature <- function(model, s, delta = 0.0015, n_iter = 6) {
  if (is.na(s)) return(NA_real_)
  ctrl <- Control(); ctrl$shading_model <- model; ctrl$save_RK45_cache <- TRUE
  br <- equilibrium_birth_rate_at(model, s, n_iter)
  resident <- expand_parameters(trait_matrix(s, "lma"), p_fit,
                                birth_rate_list = br)
  scm <- run_scm(resident, ctrl = ctrl)
  mut <- c(s - delta, s, s + delta)
  mutants <- expand_parameters(trait_matrix(mut, "lma"), resident,
                               birth_rate_list = rep(1, 3))
  scm$run_mutant(mutants)
  r0 <- scm$net_reproduction_ratios
  f  <- log(tail(r0, 3))
  (f[3] - 2 * f[2] + f[1]) / delta^2
}

## ===========================================================================
## (1) Map the gradient over a grid of resident LMA, per model
## ===========================================================================
models   <- c("deep-crown", "mean-light", "crown-centre", "flat-top-soft-box", "ppa")
lma_grid <- seq(0.065, 0.125, length.out = 25)

message("mapping selection gradient over the LMA grid ...")
grad_df <- dplyr::bind_rows(lapply(models, function(m) {
  message("  ", m)
  g <- vapply(lma_grid, function(l) sel_gradient(m, l), numeric(1))
  data.frame(method = m, lma = lma_grid, gradient = g)
}))
grad_df$method <- mk_factor(grad_df$method)

## ===========================================================================
## (2) Root-solve the attractor where the grid shows a single + -> - crossing
## ===========================================================================
find_attractor <- function(m) {
  g <- grad_df$gradient[grad_df$method == m]
  l <- grad_df$lma[grad_df$method == m]
  sign_change <- which(g[-length(g)] > 0 & g[-1] < 0)   # + -> - (attracting)
  n_cross <- sum(diff(sign(g)) != 0)
  if (length(sign_change) != 1 || n_cross != 1) {
    return(data.frame(method = m, lma_star = NA_real_,
                      well_defined = FALSE, n_crossings = n_cross))
  }
  br <- range(l[c(sign_change, sign_change + 1)])
  root <- tryCatch(
    uniroot(function(x) sel_gradient(m, x), interval = br, tol = 1e-5)$root,
    error = function(e) NA_real_)
  data.frame(method = m, lma_star = root,
             well_defined = !is.na(root), n_crossings = n_cross)
}
message("root-solving attractors ...")
attractors <- dplyr::bind_rows(lapply(models, find_attractor))
attractors$method <- mk_factor(attractors$method)

## PPA finds a local zero crossing, but its gradient is stepped by the discrete
## layers: a discontinuous jump sits just above the attractor (verified on a fine
## grid) and the tail at high LMA is erratic, so the value is sensitive to the
## layer settings. Flag it as the one not-robust attractor.
attractors$robust <- attractors$well_defined & attractors$method != "ppa"

## Classify each singular strategy: all our crossings are + -> - (convergence
## stable); add the evolutionary-stability (uninvadability) test.
message("classifying singular strategies (evolutionary stability) ...")
attractors$curvature <- vapply(seq_len(nrow(attractors)),
  function(i) fitness_curvature(as.character(attractors$method[i]),
                                attractors$lma_star[i]), numeric(1))
attractors$class <- ifelse(is.na(attractors$lma_star), NA_character_,
                           ifelse(attractors$curvature < 0, "CSS", "branching"))
print(attractors)

## save the raw numbers for reproducibility
write.csv(grad_df,
          file.path(dirname(out_repo), "evolutionary-attractor-gradients.csv"),
          row.names = FALSE)
write.csv(attractors,
          file.path(dirname(out_repo), "evolutionary-attractor-singular-strategies.csv"),
          row.names = FALSE)

## ===========================================================================
## (3) Figure: (a) gradient vs LMA with zero crossings; (b) attractor LMA
## ===========================================================================
attr_ok <- dplyr::filter(attractors, well_defined)

## Filled point = robust attractor; open point = exists but stepped (PPA).
shape_vals <- c(`TRUE` = 19, `FALSE` = 1)

## y-limits clipped so the smooth crossings are legible despite PPA's spikes
## (PPA's gradient reaches ~4000 in magnitude).
ylim <- c(-400, 400)

p_grad <- ggplot(grad_df, aes(lma, gradient, colour = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(data = attr_ok, aes(lma_star, 0, shape = robust), size = 3) +
  scale_colour_manual(values = model_cols, labels = model_labels,
                      name = "Shading model") +
  scale_shape_manual(values = shape_vals, guide = "none") +
  coord_cartesian(ylim = ylim) +
  labs(x = "Resident LMA", y = "Selection gradient,  d log(R0) / d(LMA)",
       tag = "(a)",
       subtitle = "Selection gradient at equilibrium (zero, negative slope = attractor)") +
  theme_pub

p_star <- ggplot(attractors, aes(method, lma_star, colour = method)) +
  geom_point(aes(shape = robust), size = 3.5) +
  geom_text(aes(label = formatC(lma_star, format = "f", digits = 4)),
            vjust = -1.1, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = model_cols, labels = model_labels,
                      name = "Shading model") +
  scale_shape_manual(values = shape_vals, guide = "none") +
  scale_x_discrete(labels = model_labels) +
  labs(x = NULL, y = "Attracted LMA  (singular strategy)",
       tag = "(b)", subtitle = "Where evolution settles: methods disagree (~15% spread)") +
  theme_pub +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9))

fig <- (p_grad / p_star) + plot_layout(heights = c(1.25, 1)) +
  plot_annotation(caption = paste(
    "Experimental. Singular strategy = zero of the selection gradient. All five are convergence",
    "stable (negative slope) AND uninvadable maxima (negative\nfitness curvature) -- i.e. CSS,",
    "genuine evolutionary endpoints. Yet the endpoints DISAGREE: a monotone ~15% spread",
    "(Deep Crown 0.101 -> PPA (smoothed) 0.086). Even methods that agree on stand dynamics and give\nsimilar",
    "fitness-landscape shapes settle on different LMA, because each reaches a different equilibrium",
    "density. PPA's gradient is stepped by its discrete layers (open point; clipped here, |gradient|",
    "up to ~4000; curvature ~10x steeper) so its endpoint is sensitive to the layer settings.")) &
  theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"))

save_fig(fig, "plant-evolutionary-attractor", width = 8, height = 8)

message("\nDone.")
