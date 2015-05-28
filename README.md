# The Trait Ecology and Evolution (plant) model

[![Build Status](https://travis-ci.org/traitecoevo/plant.png?branch=master)](https://travis-ci.org/traitecoevo/plant)

## Installation

Installation requires a C++11 compatible C compiler (OSX >= 10.10/Yosemite satisfies this, plus we've had success on Ubuntu 12.4 and 14.4).

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
