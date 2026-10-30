## Seat a TF24 individual at a given size and allometric departure, with reserves
## at a_st3 of ITS capacity and the tracked marginal leaf return at its own
## current value -- what a newborn gets at birth size, moved to this size.
##
## set_initial_states() alone cannot do this: it seats a NEWBORN, writing the
## seed's height over whatever height was set before it, so a plant set to 5 m
## and then seeded comes back at birth size.
##
## A copy of this lives in overstorey_staging/allometry_demo_helpers.R, which
## renders outside the test suite.
tf24_seat <- function(ind, env, height, phi = 0, psi = 0) {
  ind$set_initial_states(env)
  ind$set_state("height", height)
  ind$set_state("log_area_leaf_departure", phi)
  ind$set_state("log_area_sapwood_departure", psi)
  s <- ind$strategy
  ## mass_sapwood is linear in sapwood area, A*(h) exp(phi) theta exp(psi).
  mass_sapwood <- plant:::TF24_strategy_expand_allometry(s, height, 0, 0)$mass_sapwood *
    exp(phi) * exp(psi)
  ind$set_state("storage", s$pars$a_st3 * s$pars$a_st1 * mass_sapwood)
  ind$compute_rates(env)
  ind$set_state("leaf_marginal_return_tracked", ind$aux("leaf_marginal_return"))
  invisible(ind)
}
