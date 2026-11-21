# FF16 stands for the gradient tests, and the Control their referees run at.

# The birth-date coordinate the sweep differentiates on, at a tight seed-height
# root tolerance.
#
# ⚠️ THE REFEREE RUNS AT A TIGHT SEED-HEIGHT ROOT TOLERANCE, set here although
# it is the default, so a looser default cannot weaken the referee; the
# gradient does not need it. The seed's height is a root-find's iterate in the
# forward run and the exact root's derivative on the differentiated path
# (implicit function theorem). A difference at `offspring_height_tol = 1e-8`
# differentiates the iterate, whose stopping point moves with the trait
# differently from the root: measured on a five-year stand, d(mass)/d(lma)
# through the FF16 trade-off function was off by 6.8e-4 at a relative step of
# 1e-4 and by 1.8e-4 at 1e-5, while the sweep's answer did not move with the
# tolerance at all. At 1e-12 the two agree to 2.5e-8.
ff16_referee_control <- function(...) {
  Control(node_density_coordinate = "birth_date", offspring_height_tol = 1e-12,
          ...)
}

# FF16 at the testing parameters, one species per entry of `lma`, each at a
# constant birth rate (`birth_rate` is recycled over the species). The default
# stand reproduces: hmat is lowered to 3 m so cohorts mature inside its eight
# years, which is what puts the reproduction parameters on the gradient's path.
# `hmat = NULL` keeps the testing value, for a stand that only grows.
ff16_stand_parameters <- function(lma = 0.0825, hmat = 3, lifetime = 8,
                                  birth_rate = 2) {
  stopifnot(length(birth_rate) %in% c(1L, length(lma)))
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  m <- if (is.null(hmat)) trait_matrix(lma, "lma") else
    trait_matrix(c(lma, rep(hmat, length(lma))), c("lma", "hmat"))
  add_strategies(p, m, birth_rate = as.list(rep_len(birth_rate, length(lma))))
}

# `p` run with its states recorded, ready for stand_gradient(); by default the
# reproducing stand at the referee's Control.
ff16_stand_run <- function(p = ff16_stand_parameters(),
                           ctrl = ff16_referee_control(), ...) {
  run(p, NULL, ctrl, record = "states", ...)$solver
}

# The growing stand the census tests read: one species per `lma` at the testing
# hmat, so no cohort matures, run for `lifetime` years and recorded at the
# referee's Control.
ff16_stand <- function(lifetime = 5, lma = 0.0825, birth_rate = 2,
                       ctrl = ff16_referee_control(), ...) {
  ff16_stand_run(ff16_stand_parameters(lma, hmat = NULL, lifetime, birth_rate),
                 ctrl, ...)
}
