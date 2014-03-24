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
