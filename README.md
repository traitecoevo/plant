# plant: A package for modelling forest trait ecology and evolution

[![Build Status](https://travis-ci.org/traitecoevo/plant.png?branch=master)](https://travis-ci.org/traitecoevo/plant)
[![Build status](https://ci.appveyor.com/api/projects/status/github/traitecoevo/plant?branch=master&svg=true)](https://ci.appveyor.com/project/traitecoevo/plant/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/traitecoevo/plant/badge.svg?branch=master)](https://coveralls.io/github/traitecoevo/plant?branch=master)


The plant package for R is an extensible framework for modelling size- and trait-structured demography, ecology and evolution in simulated forests. At its core, plant is an individual-based model where plant physiology and demography are mediated by traits. Individual plants from multiple species can be grown in isolation, in patches of competing plants or in metapopulations under a disturbance regime. These dynamics can be integrated into metapopulation-level estimates of invasion fitness and vegetation structure. Accessed from R, the core routines in plant are written in C++. The package provides for alternative physiologies and for capturing trade-offs among parameters. A detailed test suite is provided to ensure correct behaviour of the code.

## Citation

Falster DS, FitzJohn RG, Brännström Å, Dieckmann U, Westoby M (2016) plant: A package for modelling forest trait ecology & evolution. *Methods in Ecology and Evolution* 7: 136-146. doi: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)

## Documentation

An overview of the plant package is given by the above publication. Further background on the default `FF16` growth model is available in Falster *et al* 2011 ([10.1111/j.1365-2745.2010.01735.x](http://doi.org/10.1111/j.1365-2745.2010.01735.x)) and Falster *et al* 2017 ([10.1101/083451](http://doi.org/10.1101/083451)).

`plant` comes with a lot of documentation, available at [https://traitecoevo.github.io/plant/](https://traitecoevo.github.io/plant/). Initial versions for some of the material there was also  included as supplementary material with the publication about plant, which can be accessed [here](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12525/abstract#footer-support-info). 


## Package structure

Plant is a complex package, using [C++11](https://en.wikipedia.org/wiki/C%2B%2B11) behind the scenes for speed with [R6 classes](https://cran.r-project.org/web/packages/R6/vignettes/Introduction.html) (via the [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html) and [RcppR6](https://github.com/richfitz/RcppR6) packages).  In this blog post, Rich FitzJohn and I describe the [key technologies used to build the plant package](https://methodsblog.wordpress.com/2016/02/23/plant/). 

If you are interested in developing plant you should read the [Developer Notes](https://traitecoevo.github.io/plant/articles/developer_notes.html).

## Installation

**Requirements**

- You must be using R 3.3.0 or newer. At this stage the package is not on CRAN. You're options for installing are described below.


- Installation requires a C++11 compatible C compiler (OSX >= 10.10/Yosemite satisfies this, as do standard linux Ubuntu 12.04 and 14.04). On Windows machines you will need to install [Rtools](http://cran.r-project.org/bin/windows/Rtools/). When I tried this in [Rstudio](https://www.rstudio.com/), the program [automagically](https://en.oxforddictionaries.com/definition/automagically) sensed the absence of a compiler and asked if I wanted to install Rtools. Click `Yes`!

- The `plant` package can be installed direct from github using the [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html) package. `plant` also requires the packages `loggr` and `RcppR6` packages. Install those with

```r
install.packages("remotes")

remotes::install_github("smbache/loggr", dependencies=TRUE)
remotes::install_github("richfitz/RcppR6", dependencies=TRUE)
```

- `plant` also depends on several packages available from CRAN. You can either install these yourself or let `remotes` handle if following installation **Option 1**. Package management in R can sometimes be messy, if you run into errors please try installing these pacakges one at a time.

```r

install.packages(c("Rcpp", "R6", "crayon", "nleqslv", "BB" ,"BH"))
``` 

**Option 1, using `remotes::install_github`**

```r
remotes::install_github("traitecoevo/plant", dependencies=TRUE)
```

or if you have already installed the `plant` dependencies

```r
remotes::install_github("traitecoevo/plant")
```

**Option 2, download and install locally**

If installing locally you will still need to install the `loggr` and `RcppR6` packages. Install using `remotes::install_github` as above or, alternatively, download a zip file from github:

- [RcppR6](https://github.com/richfitz/RcppR6/archive/master.zip)
- [loggr](https://github.com/smbache/loggr/archive/master.zip) (if needed)
- [plant](https://github.com/traitecoevo/plant/archive/master.zip). 

Unzip these archives and then for each package run the command

```r
install.packages("path_to_package", repos = NULL, type="source")

```
where `path_to_package` is the folder for each package, e.g. `~/Downloads/plant-master`


**Option 3, installing different versions** 

To install a specific (older) release, decide for the version number that you want to install in https://github.com/traitecoevo/plant/releases  e.g.

```r
remotes::install_github("traitecoevo/plant", ref = "v1.0.0", dependencies=TRUE)

```

with `"v1.0.0"` replaced by the appropriate version number. Note, the latest version of `plant` resides on the `develop` branch, which is sporadically relased. `plant` follows [semantic versioning](https://semver.org/) meaning that major version indicate a potential break in backward compatibility.

**Option 4, building from source**

If familiar with [git](https://git-scm.com/) you might find it easiest to build `plant` directly from the source code. This is most useful if developing new models or strategies, or to contribute new features.

First, clone the `plant` repository

```
git clone https://github.com/traitecoevo/plant
```

then in the terminal or command line

```
cd plant
make
```

If using Rstudio, you might like to use `devtools`
```
devtools::load_all('path/to/plant')
```

## Usage

Plant has been used in the following publications:

- Falster DS, FitzJohn RG, Brännström Å, Dieckmann U, Westoby M (2016) plant: A package for modelling forest trait ecology & evolution. *Methods in Ecology and Evolution* 7: 136-146. doi: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)&nbsp; code: [github](https://github.com/traitecoevo/plant_paper)
- Falster DS, Duursma RA, FitzJohn RG (2016) Trajectories: how functional traits influence plant growth and shade tolerance across the life-cycle. *bioRxiv*: 083451. doi: [10.1101/083451](http://doi.org/10.1101/083451)&nbsp; code: [github](https://github.com/traitecoevo/growth_trajectories)
