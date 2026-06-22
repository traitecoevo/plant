# plant: A package for modelling forest trait ecology and evolution

<!-- badges: start -->
[![R-CMD-check](https://github.com/traitecoevo/plant/workflows/R-CMD-check/badge.svg)](https://github.com/traitecoevo/plant/master)
[![Codecov test coverage](https://codecov.io/gh/traitecoevo/plant/branch/master/graph/badge.svg)](https://codecov.io/gh/traitecoevo/plant?branch=master)
<!-- badges: end -->

The plant package for R is an extensible framework for modelling size- and trait-structured demography, ecology and evolution in simulated forests. At its core, plant is an individual-based model where plant physiology and demography are mediated by traits. Individual plants from multiple species can be grown in isolation, in patches of competing plants or in metapopulations under a disturbance regime. These dynamics can be integrated into metapopulation-level estimates of invasion fitness and vegetation structure. Accessed from R, the core routines in plant are written in C++. The package provides for alternative physiology models and for capturing trade-offs among parameters. A detailed test suite is provided to ensure correct behaviour of the code.

> **Development status.** Active development happens on the `develop` branch, which is currently well ahead of the most recent tagged release. Notably, `develop` includes a new physiological model, **TF24**, which we plan to document and ship in an upcoming release. If you want the latest features (including TF24), install from `develop` (see [Installation](#installation)); for a stable, citable version, install a tagged [release](https://github.com/traitecoevo/plant/releases).
>
> Current development is tracked on our [project board](https://github.com/orgs/traitecoevo/projects/5) and in the [issue tracker](https://github.com/traitecoevo/plant/issues) — that's the best place to see what's planned, report bugs, or pick up something to contribute.

## Citation

Falster DS, FitzJohn RG, Brännström Å, Dieckmann U, Westoby M (2016) plant: A package for modelling forest trait ecology & evolution. *Methods in Ecology and Evolution* 7: 136-146. doi: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)

## Documentation

An overview of the plant package is given by the above publication. Further background on the default `FF16` growth model is available in Falster *et al* 2011 ([10.1111/j.1365-2745.2010.01735.x](http://doi.org/10.1111/j.1365-2745.2010.01735.x)) and Falster *et al* 2017 ([10.1101/083451](http://doi.org/10.1101/083451)).

The narrative documentation — user guides, theory, and worked examples — now lives at **[Overstorey](https://traitecoevo.github.io/overstorey/)**, a dedicated field guide to the `plant` model ([source](https://github.com/traitecoevo/overstorey)). The package's own [pkgdown site](https://traitecoevo.github.io/plant/) hosts the API/function reference. Initial versions of some of this material were also included as supplementary material with the publication about plant, which can be accessed [here](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12525/abstract#footer-support-info). 

## Package structure

Plant is a complex package, using [C++20](https://en.wikipedia.org/wiki/C%2B%2B20) behind the scenes for speed with [R6 classes](https://cran.r-project.org/web/packages/R6/vignettes/Introduction.html) (via the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppR6](https://github.com/richfitz/RcppR6) packages).  In this blog post, Rich FitzJohn and I describe the [key technologies used to build the plant package](https://methodsblog.wordpress.com/2016/02/23/plant/). 

If you are interested in developing or extending plant, start with [agents.md](agents.md), which documents the package architecture, the C++/R interface, the build workflow, and how to add a new model. The `new-strategy` skill (`.claude/skills/new-strategy/`) walks through scaffolding and implementing a new strategy.

## Installation

**Requirements**

- You must be using R 4.5.0 or newer. At this stage the package is not on CRAN. Your options for installing are described below.

- Installation requires a [C++20](https://en.wikipedia.org/wiki/C%2B%2B20) compatible compiler (recent versions of clang/gcc on macOS and Linux satisfy this). On Windows machines you will need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/). When I tried this in [Rstudio](https://www.rstudio.com/), the program [automagically](https://en.oxforddictionaries.com/definition/automagically) sensed the absence of a compiler and asked if I wanted to install Rtools. Click `Yes`!

**Option 1, using `remotes::install_github`**

The `plant` package can be installed direct from github using the [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html) package:

```r
remotes::install_github("traitecoevo/plant", dependencies=TRUE)
```

To install a specific (older) release, decide for the version number that you want to install in https://github.com/traitecoevo/plant/releases  e.g.

```r
remotes::install_github("traitecoevo/plant@v1.0.0", dependencies=TRUE)
```

with `"v1.0.0"` replaced by the appropriate version number. Note, the latest version of `plant` resides on the `develop` branch, which is sporadically released. `plant` follows [semantic versioning](https://semver.org/) meaning that major version indicate a potential break in backward compatibility.

**Option 2, building from source**

If familiar with [git](https://git-scm.com/) you might find it easiest to build `plant` directly from the source code. This is most useful if developing new models or strategies, or to contribute new features.

First, clone the `plant` repository

```
git clone https://github.com/traitecoevo/plant
```

Open an R session in the folder, then to install dependencies run

```
devtools::install_deps()
```

Then to compile the project

```
devtools::install()
```
or 

```
devtools::load_all()
```

## Getting started

The best place to start is **[Overstorey](https://traitecoevo.github.io/overstorey/)**, which hosts the full narrative documentation and a set of tutorial guides. In particular:

- *Get started with plant* — an overview of the package and a minimal worked example.
- *Individuals*, *Patch dynamics* and *Demography* — the core simulation workflows.
- *Parameters* — configuring strategies, traits and control settings.

The package's [pkgdown site](https://traitecoevo.github.io/plant/) hosts the function reference.

A minimal example, growing a patch of competing plants with the default `FF16` model:

```r
library(plant)

# Set up parameters for the default FF16 model
p <- scm_base_parameters("FF16")
p <- expand_parameters(trait_matrix(0.0825, "lma"), p)

# Run the deterministic (method-of-characteristics) solver and collect output
results <- run_scm_collect(p)
```

If you want to develop or extend `plant` (e.g. add a new strategy/model), see [agents.md](agents.md) and the `new-strategy` skill (`.claude/skills/new-strategy/`).

## Benchmarking

For fair performance comparisons across branches, rebuild compiled code before
running benchmarks:

```sh
make
Rscript -e "devtools::load_all(quiet=TRUE); run_plant_benchmarks()"
```

Running benchmarks without a fresh `make` can compare stale binaries and give
misleading timing differences.

## Getting help

Questions, bug reports and feature requests are welcome via the [GitHub issue tracker](https://github.com/traitecoevo/plant/issues).

## Publications using plant

Here are some example publications using plant:

- Falster DS, FitzJohn RG, Brännström Å, Dieckmann U, Westoby M (2016) plant: A package for modelling forest trait ecology & evolution. *Methods in Ecology and Evolution* 7: 136-146. DOI: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)&nbsp; code: [github](https://github.com/traitecoevo/plant_paper)
- Falster DS, Duursma RA, FitzJohn RG (2018) How functional traits influence plant growth and shade tolerance across the life cycle. *Proceedings of the National Academy of Sciences* 115: E6789–E6798. DOI: [10.1073/pnas.1714044115](http://doi.org/10.1073/pnas.1714044115)&nbsp; code: [github](https://github.com/traitecoevo/growth_trajectories)
- Falster DS, Kunstler GK, FitzJohn RG, Westoby M (2021) Emergent shapes of trait-based competition functions from resource-based models: a Gaussian is not normal in plant communities. *The American Naturalist* 198: 256–267. DOI: [10.1086/714868](http://doi.org/10.1086/714868)&nbsp; code: [github](https://github.com/traitecoevo/competition_kernels)



## Contributing

Contributions are welcome. By submitting a pull request or code to this repository, you agree to the terms of the [Contributor License Agreement](CLA.md).

