# The TRait Ecology and Evolution (TREE) model

## Installation

Installation requires a C++11 compatible C compiler (OSX >= 10.10/Yosemite satisfies this, plus we've had success on Ubuntu 12.4 and 14.4).

**Option 1, using [`drat`](https://github.com/eddelbuettel/drat)**

(install `drat` with `install.packages("drat")`)

```r
drat:::add("traitecoevo")
install.packages("tree2")
```

(versions of R before 3.2.0 may require `install.packages("tree2", type="source")` here)

**Option 2, using `devtools::install_github`**

```r
devtools::install_github("traitecoevo/tree2")
```


```
Cohort:
  name_cpp: "tree2::Cohort<T>"
  templates:
    parameters: T
    concrete:
      - ["FFW16": "tree2::FFW16_Plant"]
  constructor:
    name_cpp: "tree2::make_cohort<T>"
    args: [strategy: "T::strategy_type"]
  active:
    plant: {type: "T", access: field, readonly: true}
    height: {type: double, access: member}
    area_leaf: {type: double, access: member}
    fecundity: {type: double, access: member}
    ode_size: {type: size_t, access: member}
    ode_state: {type: "tree2::ode::state_type", access: function, name_cpp: "tree2::ode::r_ode_state", name_cpp_set: "tree2::ode::r_set_ode_state"}
    ode_rates: {type: "tree2::ode::state_type", access: function, name_cpp: "tree2::ode::r_ode_rates"}
  methods:
    area_leaf_above:
      return_type: double
      args: [height: double]
    growth_rate_gradient:
      return_type: double
      name_cpp: "r_growth_rate_gradient"
      args: [environment: "const tree2::Environment&"]
    compute_vars_phys:
      return_type: void
      args: [environment: "const tree2::Environment&"]
    compute_initial_conditions:
      return_type: void
      args: [environment: "const tree2::Environment&"]
```
