# The Trait Ecology and Evolution (plant) model

[![Build Status](https://travis-ci.org/traitecoevo/plant.png?branch=master)](https://travis-ci.org/traitecoevo/plant)

## Documentation

`plant` comes with a lot of documentation.  Run `library(help=plant)` to see the index, and see vignettes:

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

## Installation

Installation requires a C++11 compatible C compiler (OSX >= 10.10/Yosemite satisfies this, plus we've had success on Ubuntu 12.04 and 14.04).

**Option 1, using [`drat`](https://github.com/eddelbuettel/drat)**

(install `drat` with `install.packages("drat")`)

```r
drat:::add("traitecoevo")
install.packages("plant", dependencies=TRUE)
```

Versions of R before 3.2.0 may require `install.packages("plant", type="source")` here.  The `dependencies=TRUE` installs some additional packages that are used in tests and some of the routines.

**Option 2, using `devtools::install_github`**

```r
devtools::install_github("traitecoevo/plant", dependencies=TRUE)
```

Development versions of `plant` (version 0.2.3 onwards) depends on the `loggr` package: install that from drat or with

```
devtools::install_github("smbache/loggr")
```

### Windows

Installation on Windows is likely to be a challenge, because of the lack of a C++11 compiler.  The current Windows [toolchain](http://cran.r-project.org/bin/windows/Rtools/) uses gcc 4.6.3, which des not feature enough C++11 support to successfully compile plant.  There is support coming for gcc 4.9.2, but that is [under development](https://rawgit.com/kevinushey/RToolsToolchainUpdate/master/mingwnotes.html).  This should become available on Windows very soon (it was initially aimed for 3.2.0 but didn't quite make it).  As soon as this is available we will provide Windows binaries.

## Building vignettes

See [`docs/README.md`](docs/README.md) for details.  The vignettes are not built as part of the package installation as they take a couple of hours to build.
