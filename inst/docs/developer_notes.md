---
title: "Developer Notes"
---

## Issues

We use [github issues](https://github.com/traitecoevo/plant/issues/) for feature requests, bug reports and  developer coordination.

## Gitflow workflow

We use the [Gitflow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) workflow.
This means we have two main branches:

1. a master branch (`master`) for stable releases
2. a development branch (`develop`) where all development takes place

We also use feature branches which branch from the development branch where we can then introduce new features and then pull them back into the development branch.

##  Code Structure
#### Directory  Structure
```sh
├── R # R functions
├── docker # Docker files
│   └── plant-test
├── inst
│   ├── docs # Files for building the documentation
│   │   ├── R
│   │   ├── figure
│   │   ├── figures
│   │   ├── reference
│   │   ├── src
│   ├── include
│   │   ├── plant # .h files for plant.
│   │   └── tk
│   └── reference_plant_ff16 # Reference plant physiology ff16.
│       ├── R
│       └── falster-traitdiversity
├── man # R Documentation ie .Rd files.
├── scripts # R scripts
├── src # C++ files
├── tests # Tests
└── vignettes # Empty until the documentation is built
```

## Building the website:

To build and push the website to the `gh-pages` branch run `make website`.
This first builds the vignettes from `/inst/docs` then builds the site with [pkgdown](http://pkgdown.r-lib.org).
The site is then pushed to the `gh-pages` branch. (https://traitecoevo.github.io/plant)

## Strategies

The foundation of the `plant` package are sub-models for an individual species’ physiological strategy. The `FF16` strategy is the default strategy. However, the plant package has been written to enable new physiological strategy to be added, and then use the same machinery for modelling size-structured population dynamics. This possible via the use of templating in C++ 11. 

Adding a new physiological strategy requires some changes to be made throughout the code and then the code recompiled. To achieve this, you'll need to download the package code from Github (new strategies cannot be added using an installed R package) and follow the instructions below. 

As a naming convention, strategies are named using two initials r (eg. the first two authors) followed by a year, e.g. `FF16`. A single letter suffix can be added to indicate a minor modification of an existing strategy, e.g. `FF16r`. 

### Installing new strategies

Let's assume you want to add a new strategy called `XX99`. To add this strategy, the following changes are required. First, some new files need to be created:

* `inst/include/plant/XX99_strategy.h` (C++ header file, describing strategy)
* `src/XX99_strategy.cpp`  (C++ source file, describing strategy)
* `R/XX99.R`  (R wrapper functions for the new strategy, plus hyperparameter function)
* `tests/testthat/test-strategy-XX99.R` (Tests of the new strategy)

In addition, the following files must be modified to include appropriate details for the new strategy (use the examples for `FF16` as a guide):

* `inst/include/RccpR6_classes.yml`
* `inst/include/plant.h`
* `src/plant_plus.cpp`
* `src/plant_tools.cpp`
* `tests/testthat/helper-plant.R`
* `R/scm_support.R` (Currently uses the FF16 hyperpar and make_hyperpar)

Rather than making the above modifications by hand, you can use the scaffolder found in `inst/scripts/new_strategy_scaffolder.R` to create a suitable template. The scaffolder is set up to add code at appropriate points. 

Before trying to install a new strategy, make sure your project is at a state where you are happy to go back to if you decide to, ie: do this in a new branch or fork. There is no easy way to undo the changes so you will have to do `git reset --hard HEAD` if you want to go back. You will still have to delete the new and untracked files listed above.

To install a new strategy using the scaffolder, run the following code from an R session at the project root:

```{r}
source('inst/scripts/new_strategy_scaffolder.R')
create_strategy_scaffold("XX99", "FF16")
```
where `XX99` is the name of your new strategy and `FF16` is the strategy that you want to copy from. If run without `FF16`, the scaffolder copies from the `FF16` strategy by default.

This will create the code-base for a new template. Then run `make clean; make; make test`to recompile and run the test suite with your new strategy template.

If that works you can then modify the files `src/XX99_strategy.cpp` and `inst/include/plant/XX99_strategy.h` to reflect any changes in biology you desire for the new strategy. If you make changes to your new strategies parameters you will have to update the files: 
	* `inst/include/RccpR6_classes.yml` - see top-level section under the new strategy name, describing your parameters.
	* `tests/testthat/test-strategy-XX99.R` with suitable tests for your new strategy.

Then recompile and test as above.

After the project compiles, you will notices even more files that have changed, for example in files like `NAMESPACE`, `R/RcppExports.R` etc. These files are auto-generated and therefore modifications are triggered from files listed above.

----

## Makefile guide:

#### All

```make all```

Compiles the package and builds the function documentation to `/man`. 

#### Test

```make test```

Runs the test suite.

#### Install

```make install```

Installs this package (` R CMD INSTALL .`).

#### Build

``` make build```

Builds plant (` R CMD build --no-build-vignettes`).

#### Check

```make check```

[Checks](http://r-pkgs.had.co.nz/check.html) plant for common problems.

#### Clean

```make clean```

Deletes `*.o` and `*.so` files that were produced when compiling the package.

#### vignettes

```make vignettes```

Builds the vignettes.

#### pkgdown

```make pkgdown```

Makes the pkgdown site.

#### website

```make website```

Pushes the site to the `gh-pages` branch. 
