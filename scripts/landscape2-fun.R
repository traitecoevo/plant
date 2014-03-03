# This bit is way uglier than needed; need to expand out a schedule to
# include mutants.
expand.schedule <- function(schedule, n.mutant) {
  n.resident <- schedule$n_species
  ret <- new(CohortSchedule, n.resident + n.mutant)
  ret$max_time <- schedule$max_time

  ## Copy residents over:
  for (i in seq_len(n.resident))
    ret$set_times(schedule$times(i), i)

  ## Introduce mutants at all unique times:
  times.mutant <- unique(sort(unlist(schedule$all_times)))
  for (i in seq_len(n.mutant))
    ret$set_times(times.mutant, n.resident + i)

  ret$ode_times <- schedule$ode_times

  ret
}

seq.log <- function(from, to, length.out)
  exp(seq(log(from), log(to), length.out=length.out))

landscape <- function(lma, p, schedule) {
  p.with.mutants <- p$copy()
  for (i in lma)
    p.with.mutants$add_strategy_mutant(new(Strategy, list(lma=i)))
  schedule.with.mutants <- expand.schedule(schedule, length(lma))
  ebt.with.mutants <- run.ebt(p.with.mutants, schedule.with.mutants)
  w.with.mutants <- ebt.with.mutants$fitnesses
  w.with.mutants[-seq_len(p$size)]
}

landscape.empty <- function(lma, p, schedule) {
  p.empty <- p$copy()
  p.empty$clear()
  for (i in lma)
    p.empty$add_strategy_mutant(new(Strategy, list(lma=i)))

  schedule.empty <- new(CohortSchedule, length(lma))
  schedule.empty$max_time  <- schedule$max_time
  schedule.empty$ode_times <- schedule$ode_times
  schedule.empty$all_times <-
    rep(list(unique(sort(unlist(schedule$all_times)))), length(lma))
  ebt.empty <- run.ebt(p.empty, schedule.empty)
  ebt.empty$fitnesses
}
