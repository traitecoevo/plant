// -*-c++-*-
#ifndef PLANT_PLANT_RUNNER_H_
#define PLANT_PLANT_RUNNER_H_

#include <plant/node_schedule.h>

namespace plant {

// Shared skeleton for the two schedule-driven patch runners:
//   * SCM<T,E>                 -- the deterministic (method-of-characteristics)
//                                 solver, which introduces cohorts (Nodes); and
//   * StochasticPatchRunner<T,E> -- the finite-population solver, which
//                                 introduces individuals and kills them as
//                                 discrete events.
//
// Both own a NodeSchedule and advance a patch by repeatedly consuming the next
// scheduled event until the schedule is exhausted. That common lifecycle lives
// here so the two runners present a *consistent interface* and don't duplicate
// the schedule plumbing, the completion test, and the run loop.
//
// This is a CRTP base (static, compile-time polymorphism) rather than a virtual
// base on purpose: the runners sit just above the model's hot path and the
// codebase deliberately avoids runtime dispatch here (see agents.md s12). The
// per-event work that genuinely differs between the two solvers -- batched
// cohort introduction plus multiple integration modes for the SCM, versus a
// single arrival plus a stochastic death step for the stochastic runner -- is
// left to the derived class via the contract below, so it stays inlinable.
//
// Contract the derived class must provide (called through the CRTP downcast):
//   void  reset();      -- return patch / schedule / solver to their t=0 state
//   <any> run_next();   -- consume the next scheduled event (return ignored)
//
// The base owns the schedule (named `node_schedule` to match the field both
// solvers historically used) and implements complete()/run() in terms of it.
// run() is a sensible default (reset, then step to completion); a derived class
// that needs extra per-step bookkeeping -- as the SCM does, to collect history
// and refinement errors -- simply defines its own run() to hide this one.
template <typename Derived>
class ScheduleDrivenRunner {
public:
  // Run the whole schedule from a fresh reset to completion.
  void run() {
    Derived& self = derived();
    self.reset();
    while (!complete()) {
      self.run_next();
    }
  }

  // True once every scheduled event has been consumed.
  bool complete() const { return node_schedule.remaining() == 0; }

protected:
  // The derived class builds the schedule from its Parameters and hands it in.
  explicit ScheduleDrivenRunner(NodeSchedule s) : node_schedule(s) {}

  Derived&       derived()       { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }

  NodeSchedule node_schedule;
};

} // namespace plant

#endif
