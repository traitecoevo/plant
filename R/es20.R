##' Create an ES20 Individual
##'
##' ES20 is FF16 with a non-structural carbohydrate store (Stefaniak et al.,
##' egusphere-2026-1474). Unlike the other models, an ES20 plant carries its
##' tissues and its storage pool as states, so a new individual is seeded from
##' the seedling the strategy describes (\code{set_initial_states}); nodes and
##' stochastic arrivals in a patch are seeded the same way by the core.
##' @title Create an ES20 Individual
##' @param s An \code{\link{ES20_Model}} object
##' @param env The environment the seedling is born into
##'   (\code{\link{ES20_Environment}})
##' @export
##' @rdname ES20_Individual
##' @examples
##' pl <- ES20_Individual()
##' pl$state("mass_storage")
ES20_Individual <- function(s = ES20_Model(), env = ES20_Environment()) {
  ind <- Individual("ES20", "FF16_Env")(s)
  ind$set_initial_states(env)
  ind
}

##' @param ... Passed to the \code{Parameters} constructor (see
##'   \code{\link{FF16_Parameters}}).
##' @export
##' @rdname ES20_Individual
ES20_Parameters <- function(...) {
  Parameters("ES20", "FF16_Env")(...)
}

##' An environment for ES20: FF16's light environment, with the annual stress
##' schedule as the \code{"stress_onset"} extrinsic driver.
##'
##' Each year has a stress period running from that year's onset to the end of
##' the year, during which ES20 plants do not photosynthesise or establish and
##' respire at \code{r_min} of their usual rate. \code{stress_onset[i]} is the
##' onset in year \code{i - 1} (as a fraction of the year, so the stress lasts
##' \code{1 - stress_onset[i]} years). It is stored as a driver with one knot per
##' year at whole years, and read at the start of each year. Supply at least as
##' many years as the run spans. With no driver (\code{stress_onset = NULL})
##' there is no stress at all, which differs slightly from an onset of 1: the
##' stress factor's logistic transition still dips at the year's end.
##'
##' @title Environment for the ES20 model
##' @param stress_onset Stress onset in each year, or \code{NULL} for none.
##' @return An \code{FF16_Environment}.
##' @seealso \code{\link{ES20_stress_onset}}
##' @export
##' @examples
##' env <- ES20_Environment(ES20_stress_onset(100, duration_mean = 0.25,
##'                                           duration_sd = 30 / 365))
ES20_Environment <- function(stress_onset = NULL) {
  env <- FF16_Environment()
  if (!is.null(stress_onset)) {
    if (!is.numeric(stress_onset) || length(stress_onset) < 2L ||
        anyNA(stress_onset)) {
      stop("`stress_onset` must be a numeric vector of two or more years, ",
           "with no NA")
    }
    env$extrinsic_drivers_set_variable("stress_onset",
                                       seq_along(stress_onset) - 1,
                                       stress_onset)
  }
  env
}

##' Draw a sequence of annual stress onsets for ES20
##'
##' Each year's stress duration is drawn from a normal distribution, and the
##' onset is one minus the duration, as in Stefaniak et al.'s treatments (mean
##' durations of 54.75 and 91.25 days, standard deviations of 0, 15, 30 and 60
##' days). Draws are not truncated, as they were not in the paper.
##' @title Annual stress onsets for ES20
##' @param n_years Number of years.
##' @param duration_mean,duration_sd Mean and standard deviation of the stress
##'   duration, in years.
##' @return A numeric vector of \code{n_years} onsets, each a fraction of the
##'   year.
##' @export
ES20_stress_onset <- function(n_years, duration_mean, duration_sd = 0) {
  1 - stats::rnorm(n_years, duration_mean, duration_sd)
}

##' @param results A tidied run (see \code{\link{expand_state}}).
##' @export
##' @rdname ES20_Individual
ES20_expand_state <- function(results) {
  ## ES20 carries every tissue, the stem area and diameter as states, so there
  ## is nothing to derive from height.
  results
}
