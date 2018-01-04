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
