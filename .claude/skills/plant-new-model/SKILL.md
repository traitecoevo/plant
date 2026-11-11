---
name: plant-new-model
description: >-
  Add a new strategy model to the plant C++/R package, usually as a subclass of
  an existing strategy (as TF24f, TF24t and TF24w are of TF24), and wire it
  through the build. Use when asked to add a new model or strategy, a variant of
  an existing model, or a strategy with its own environment. Lists every file a
  new model touches. Does NOT cover changing the numerics of an existing model.
---

# Adding a strategy model to plant

A **strategy** is plant's biological model for one species (growth, mortality, reproduction, competition). An **environment** is the resource field a strategy reads (light, and for TF24 the soil column). The C++ core is templated on the `<Model, Environment>` pair, so a new model has to be named in every templated binding and every R dispatch table. There is no generator for this: the checklist in Step 2 is the wiring. The commit that added TF24t is a worked example of it (`git show --stat 0b68e62f`, then `git show 0b68e62f -- <file>`), with one caveat: it predates #731, which moved trade-offs onto `Parameters`, so its `*_hyperpar` exports and the `make_hyperpar`/`hyperpar` arms in `R/strategy_support.R` no longer exist and must not be copied. Its trade-off wiring today is the `R/default_parameters.R` row below.

Read [agents.md](../../../agents.md) §7 first. Never hand-edit generated files (`R/RcppR6.R`, `R/RcppExports.R`, `src/RcppR6.cpp`, `src/RcppExports.cpp`, `inst/include/plant/RcppR6_pre.hpp`, `RcppR6_post.hpp`, `NAMESPACE`, `man/`).

## Which route

| Route | When | In the tree |
|---|---|---|
| **Subclass** (the usual one) | A variant that keeps most of an existing strategy's biology and competes in its environment. | `TF24f`, `TF24t`, `TF24w`, each a subclass of `TF24_Model` in `TF24_Environment`. |
| **Own strategy** | New biology throughout. | `FF16`, `K93`, `TF24`. `K93` is the smallest complete model to read. |
| **Own environment** | The model reads a different resource field. | `FF16_Environment`, `K93_Environment`, `TF24_Environment`. |

Naming: two initials and a year (`FF16`), and a single-letter suffix for a variant (`TF24t`).

## Step 1 — the strategy header

For a subclass, write `inst/include/plant/models/<name>_model.h` along the lines of [tf24t_model.h](../../../inst/include/plant/models/tf24t_model.h):

```cpp
template <typename S = double>
class XX24v_Model : public TF24_Model<S> {
public:
  typedef TF24_Model<S> base_type;
  typedef std::shared_ptr<XX24v_Model> ptr;
  XX24v_Model();                       // sets this->name = "XX24v"
  static constexpr int variant_revision = 1;
};
template <typename S>
typename XX24v_Model<S>::ptr make_strategy_ptr(XX24v_Model<S> s);
```

What makes it work, and what breaks it:

- **Reuse the parent's `_Pars`.** The subclass's yml block lists `- pars: "plant::TF24_Pars<double>"` and no `<Name>_Pars` of its own exists. This is required: inherited methods read the inherited `pars`, so a shadowing struct would feed them the wrong values. New parameters go at the top level of the strategy, as TF24t's thermal traits and TF24f's acclimation settings do.
- **Overrides resolve on the concrete type.** `Individual<Sub, Env>` holds a `Sub::ptr`, so inherited methods are reused, overrides win, and the statics `state_size()` and `state_names()` resolve on the subclass. A subclass that adds states appends them after the parent's (`base_type::state_size() + 4` in TF24t) and checks the layout in `refresh_indices()`.
- **`virtual` only where a reused base method calls the step you replace.** `TF24_Model::net_mass_production_dt`, `solve_leaf` (which TF24f overrides), `fraction_allocation_reproduction`, `soil_water_potential_for_leaf` and `aux_names` are virtual for that reason. Making a base method virtual is a change to the parent and must leave its numbers bit-identical.
- **`set_initial_states()` is not virtual.** A subclass that defines it calls `base_type::set_initial_states()` first, or the parent's states (the storage pool among them) are never seeded, as TF24f, TF24t and TF24w do.
- **A subclass that is not differentiable deletes `rebind()`** (`template <class U> void rebind() const = delete;`, as TF24t does). The inherited `TF24_Model::rebind<U>()` returns the parent's type, so a gradient entry point would otherwise silently run the parent model; deleting it makes that a compile error.
- **A compound scientific version.** A subclass reports `"<parent version>.<revision>"` (TF24t is `TF24t@v<TF24 version>.<thermal_revision>`), so a change to the parent invalidates it too. See agents.md §7.

For an own strategy or environment, copy the nearest model's header and its `_Pars` struct and rename throughout; the Model contract in agents.md §2 lists what the core calls. An own strategy also needs its own `<Name>_Pars:` block in `inst/RcppR6_classes.yml`, and an own environment its own block and `templates:` entries.

## Step 2 — wire it through

Every file below names each model by hand. Search each one for the parent's name (`TF24f` is a good anchor for a TF24 subclass) and add the new model beside it.

| File | What to add |
|---|---|
| `inst/include/plant.h` | `#include <plant/models/<name>_model.h>` |
| `inst/RcppR6_classes.yml` | The `<Name>_Model:` block, and the new pair in every `templates:` list that names the parent (eleven for TF24t: `IndividualRunner`, `Individual`, `Parameters`, `Node`, `Species`, `Patch`, `SCM`, `StochasticSpecies`, `StochasticPatch`, `StochasticPatchRunner`, `OdeRunner`). An own environment adds its own block. |
| `src/<name>_node.cpp` | The `schedule_default` and `make_schedule` exports for `Parameters<<Name>, <Env>>` (copy `src/tf24t_node.cpp`). |
| `src/individual_runner.cpp` | `<Name>_oderunner_individual_internals`. |
| `src/strategy_expand.cpp` | `<Name>_strategy_expand_allometry`. |
| `src/strategy_version.cpp` | The model's arm in `strategy_scientific_version()`, which R's `model_version()` reads. |
| `R/<name>.R` | `<Name>_Individual()`, `<Name>_Parameters()` and `<Name>_expand_state()` (copy `R/tf24t.R`). |
| `R/strategy_support.R` | Every `switch()` that has an arm per model: `environment_type()`, `Environment()`, `expand_state`. The node-schedule helpers need no arm: they are found by the names `src/<name>_node.cpp` exports. |
| `R/default_parameters.R` | The `default_hyperpar()` arm, and for an own strategy its entry in the `attr(f, "traits")` switch below it, whose fallback lists TF24's traits. A TF24 subclass shares `frozen_TF24_hyperpar`. A new trade-off function belongs in the fieldtrial package, with only a frozen testing copy here. |
| `R/stand_gradient.R` | The census arm: declare one, or refuse the model with a message saying why. |
| `tests/testthat/test-model-version.R` | The model in `models`, a test of its compound version, and any top-level scientific fields in `top_level_science`; then accept the new snapshot in `_snaps/model-version.md`. |
| `tests/testthat/test-model-<name>.R` | The model's own tests. |
| `tests/testthat/test-trait-names.R`, `test-default-parameters.R` | The model in their per-model loops. |
| `tests/testthat/helper-plant.R` | `get_list_of_strategy_types()` and `get_list_of_environment_types()`, which the generic strategy tests loop over. A subclass whose appended states need seeding is left out, as TF24f and TF24t are, and tested in its own file. |

## Step 3 — regenerate, build, test

```sh
make RcppR6         # regenerates the bindings from the yml and re-adds plant's guards
make roxygen        # NAMESPACE and man/
R CMD INSTALL --no-docs -l <scratch-lib> .
R_LIBS=<scratch-lib> make test-installed
```

Use `make RcppR6` and not RcppR6 directly: its script re-adds the guards a bare regeneration drops. Then run a stand end to end:

```r
p <- default_parameters("XX24v", quiet = TRUE) |>
  add_strategies(trait_matrix(0.0825, "lma"), birth_rate = 20)
res <- run_scm(p, collect = TRUE, refine_schedule = TRUE)
```

## Verification

- A subclass that overrides nothing yet reproduces its parent's run to the bit. Check this before writing any biology; it proves the wiring.
- No existing model's number moves: `default_baseline.rds`, `reference-gradient.tsv`, `FF16_reference/` and the model-version snapshot of every other model are unchanged.
- `model_id("<Name>")` returns the compound version you intended.
- The model is refused, with a message, everywhere it is not supported (the gradient sweep, the stochastic patch, the continuous patch); TF24w's refusals are listed in agents.md §7 and TF24t's in `.claude/hazards.md`, hazard 16.
