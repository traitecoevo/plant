# Reverse-pass cost (#677)

What a `stand_gradient()` sweep costs against the run it differentiates, what has been done about it, and what is left. The issue holds the measurements; this note holds what constrains later work.

## Where it stands

odelia 0.11.0 and this branch: on the #677 FF16 stand (lifetime 60, birth-date coordinate, one species, every column) a one-metric gradient takes 0.142 s and a five-metric one 0.233 s, against a 0.013 s run. TF24 at lifetime 6 with the full schedule takes 1.54 s and 2.01 s against a 0.79 s run. Forward runs are bit-identical to odelia 0.10.5; gradients agree with it to 1e-13 relative.

## Invariants a later edit could break

- **`Patch::assign_from(from)` must leave the patch equal to `from.rebind<value_type>()` in every member a recording reads, for a `from` built from the same parameters (the patch it was rebound from).** It does not rebuild the parameters, node schedule or disturbance regime. odelia's sweep holds one active patch and re-seats it at every width through this hook (`AssignsFrom` in odelia's `ode_interface.hpp`). `rebind()` constructs a fresh patch and calls `assign_from`, so a member added to one is added to both. A member that a rate evaluation writes and a rebind would reset, but `assign_from` leaves alone, gives finite wrong gradients and raises nothing. To check a change, build once with the hook renamed, so odelia falls back to a rebind per width, and compare `stand_gradient()` bit for bit.
- **Each species' strategy is re-assigned in `assign_from`, not kept.** A TF24 strategy holds a `phylloptim::Leaf` reused across evaluations; a rebind copies the passive strategy's, and so does `assign_from`.
- **The environment is rebound afresh in `assign_from`**, so no cache a previous width's evaluations left (TF24's soil-potential and driver memos, the clamp counters) reaches the next.
- **The light spline's span coefficients are double at every scalar** (odelia's `hermite_spline`). A read at an active scalar is the passive read with the knot data's partials attached. Hand the spline its knot vectors by `std::move`: a copy is a statement per knot.
- **`XAD_REDUCED_MEMORY` is set in odelia's vendored `Config.hpp`.** plant compiles odelia's tape layout into its own objects, so plant must build against the odelia it links.

## Decided against, and why

- **Moving stage states into the System** (spike: about 1%). XAD's move leaves the value and drops the slot, so a System that reads an iterator entry twice gets a value with no derivative and nothing reports it.
- **Field-reduction rows** (spike: 2% at one metric, 6% at five). They are a second copy of `Species::field_splits_nodes` and of `CanopyShape::Q_and_q` that must agree with the first bit for bit.
- **Building the excluded field only below the boundary node** (spike: 1%). It threads a parameter through `Patch`, three environments and `ResourceSpline`.
- **A Gram–Schmidt basis of output adjoints.** The only dependence seen is one row proportional to another, which `adjoint_multiples` catches in a pairwise test.

## Open, roughly by size

1. **The light-field build at each stage.** It is still recorded operation by operation: two builds per stage (with and without the boundary interval) and a boundary-node evaluation in each. A hand-written adjoint for the build, as the crown integral has, is the largest lever left. It must not duplicate `field_splits_nodes`' arithmetic; see the decision above.
2. **Recording on worker threads.** A step's recording reads only the recorded forward state, so steps can be recorded ahead of the sweep: 2.8-3.1x on 4 threads in the spike. Needs thread-local XAD tapes (plant and odelia build with `XAD_NO_THREADLOCAL`), no function-local static scratch, and errors raised again on the main thread.
3. **The Runge-Kutta stage combinations are recorded**, though they are fixed linear algebra: 5-8% in the spike.
4. **Registering only the traits asked for:** 3.5-5.6% in the spike, exact on the kept columns.
