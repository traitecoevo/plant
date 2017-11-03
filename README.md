# plant: A package for modelling forest trait ecology and evolution

[![Build Status](https://travis-ci.org/traitecoevo/plant.png?branch=master)](https://travis-ci.org/traitecoevo/plant)
[![Build status](https://ci.appveyor.com/api/projects/status/github/traitecoevo/plant?branch=master&svg=true)](https://ci.appveyor.com/project/traitecoevo/plant/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/traitecoevo/plant/badge.svg?branch=master)](https://coveralls.io/github/traitecoevo/plant?branch=master)


The plant package for R is an extensible framework for modelling size- and trait-structured demography, ecology and evolution in simulated forests. At its core, plant is an individual-based model where plant physiology and demography are mediated by traits. Individual plants from multiple species can be grown in isolation, in patches of competing plants or in metapopulations under a disturbance regime. These dynamics can be integrated into metapopulation-level estimates of invasion fitness and vegetation structure. Accessed from R, the core routines in plant are written in C++. The package provides for alternative physiologies and for capturing trade-offs among parameters. A detailed test suite is provided to ensure correct behaviour of the code.

## Citation

Falster DS, FitzJohn RG, Brännström Å, Dieckmann U, Westoby M (2016) plant: A package for modelling forest trait ecology & evolution. *Methods in Ecology and Evolution* 7: 136-146. doi: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)

## Documentation

For an overview of the plant package see the above publication. Further background on the default `FF16` growth model is available in Falster *et al* 2011 ([10.1111/j.1365-2745.2010.01735.x](http://doi.org/10.1111/j.1365-2745.2010.01735.x)) and Falster *et al* 2017 ([10.1101/083451](http://doi.org/10.1101/083451)).

Run `library(help=plant)` to see the index of functions.

`plant` comes with a lot of documentation. Initial versions were included as supplementary material with the publication about plant, which can be accessed [here](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12525/abstract#footer-support-info). Updated versions can be built from within the package (instructions below) or accessed online:

**Details of the modelling approaches:**

* `vignette("demography")`: Modelling demography of plants, patches and metapopulations [online](https://traitecoevo.github.io/plant/vignettes/demography.pdf)
* `vignette("physiology")`: Plant physiological model [online](https://traitecoevo.github.io/plant/vignettes/physiology.pdf)

**Details of using `plant` from R:**

* `vignette("plant")`: Plant level properties [online](https://traitecoevo.github.io/plant/vignettes/plant.html)
* `vignette("cohort_spacing")`: The cohort spacing algorithm [online](https://traitecoevo.github.io/plant/vignettes/cohort_spacing.html)
* `vignette("equilibrium")`: Finding demographic equilibrium [online](https://traitecoevo.github.io/plant/vignettes/equilibrium.html)
* `vignette("patch")`: Patch level dynamics [online](https://traitecoevo.github.io/plant/vignettes/patch.html)
* `vignette("emergent")`: Patch level emergent properties [online](https://traitecoevo.github.io/plant/vignettes/emergent.html)
* `vignette("fitness")`: Calculating fitness [online](https://traitecoevo.github.io/plant/vignettes/fitness.html)
* `vignette("parameters")`: Modifying parameters of the physiological model [online](https://traitecoevo.github.io/plant/vignettes/parameters.html)

If you want to build the vignettes locally, see [`docs/README.md`](docs/README.md) for details.  The vignettes are not built as part of the package installation as they take a couple of hours to build.

**Details on package structure**

Plant is a complex package. It uses [C++11](https://en.wikipedia.org/wiki/C%2B%2B11) behind the scenes for speed, and [R6 classes](https://cran.r-project.org/web/packages/R6/vignettes/Introduction.html). 

In this blog post, Rich FitzJohn and I describe the [key technologies used to build the plant package](https://methodsblog.wordpress.com/2016/02/23/plant/).


## Installation

You must be using R 3.2.0 or newer. At this stage the package is not on cran. You're options for installing are described below.

Installation requires a C++11 compatible C compiler (OSX >= 10.10/Yosemite satisfies this, as do standard linux Ubuntu 12.04 and 14.04). On Windows machines you will need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/). When I tried this in Rstudio, Rstudio automagically sensed the absence of a compiler and asked if I wanted to install Rtools. Click `Yes`!

**Option 1, using `devtools::install_github`**

(install `devtools` with `install.packages("devtools")`)

The `plant` package can be installed direct from github using the [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html) package. `plant` also requires the packages `loggr` and `RcppR6` packages. Install those with

```r
devtools::install_github("smbache/loggr", dependencies=TRUE)
devtools::install_github("richfitz/RcppR6", dependencies=TRUE)
```

Then install plant:

```r
devtools::install_github("traitecoevo/plant", dependencies=TRUE)
```

**Option 2, download and install locally**

If installing locally you will still need to install the `loggr` and `RcppR6` packages. Install using `devtools::install_github` as above, or alternatively do as follows.

Download a zip file from github of the dependencies [RcppR6](https://github.com/richfitz/RcppR6/archive/master.zip) and [loggr](https://github.com/smbache/loggr/archive/master.zip) (if needed), plus [plant](https://github.com/traitecoevo/plant/archive/master.zip). 

Unzip these archives and then for each package run the command

```r
install.packages("path_to_package", repos = NULL, type="source")

```
where `path_to_package` is the folder for each package, e.g. `~/Downloads/plant-master`

## Usage

Plant has been used in the following publications:

- Falster DS, FitzJohn RG, Brännström Å, Dieckmann U, Westoby M (2016) plant: A package for modelling forest trait ecology & evolution. *Methods in Ecology and Evolution* 7: 136-146. doi: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)&nbsp; code: [github](https://github.com/traitecoevo/plant_paper)
- Falster DS, Duursma RA, FitzJohn RG (2016) Trajectories: how functional traits influence plant growth and shade tolerance across the life-cycle. *bioRxiv*: 083451. doi: [10.1101/083451](http://doi.org/10.1101/083451)&nbsp; code: [github](https://github.com/traitecoevo/growth_trajectories)
