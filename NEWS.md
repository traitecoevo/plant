## Plant 1.2.0  Release Notes

v1.2.0 was released on 20/03/2018

### Major Changes

- New strategy scaffolder ensures higher level operations work across all strategies and for new strategies. This includes `lcp_whole_plant` `XXX_PlantPlus`, `grow_plant_to_size`, `grow_plant_to_height`, `grow_plant_to_time`
- website now builds via `pkgdown` package (used to use staticdocs but this is deprecated)
- simplified workflow for building website
- converted supporting materials from tex to Rmd

### Minor Changes

- start documenting notes for developers in `inst/docs/developer_notes.Rmd`
- add `CITATION` file
- Address many issues in documentation and package setup causing rcmdcheck to fail

A full account of changes from the previous version is available on Github: [v1.1.0...v1.2.0](https://github.com/traitecoevo/plant/compare/v1.1.0...v1.2.0)

## Plant 1.1.0  Release Notes

v1.1.0 was released on 2/02/2018

### Major Changes

- Now compiles and runs on Windows machines (requires R 3.3.0 or newer)
- Further details on installation 
- Enable assembly_parameters to accept more named arguments 

### Minor Changes

- Add Appveyor for build tests on Windows machines
- Update tests to use latest version of testthat
- Remove package traitecoevo/callr, previously used to make system calls  
- Makefile: Add roxygen and RcppR6 to compile target 
- roxygen & Rcpp updates
- Added a `NEWS.md` file to track changes to the package.

A full account of changes from the previous version is available on Github: [v1.0.0...v1.1.0](https://github.com/traitecoevo/plant/compare/v1.0.0...v1.1.0)

_______________________________________________________________________________

## Plant 1.0.0 Release Notes

v1.0.0 was released on 23/02/2016
 
This version corresponds to the paper describing the package:

Falster, DS, RG FitzJohn, Å Brännström, U Dieckmann, M Westoby (2016) plant: A package for modelling forest trait ecology and evolution. Methods in Ecology and Evolution 7: 136-146, doi: [10.1111/2041-210X.12525](http://doi.org/10.1111/2041-210X.12525)

A full account of changes from the previous version is available on Github: [v0.2.2...v1.0.0](https://github.com/traitecoevo/plant/compare/v0.2.2...v1.0.0)

_______________________________________________________________________________


## Plant 0.2.2 Release Notes

v0.2.2 was released on 1/06/2015

Draft paper about package submitted to Methods in Ecology & Evolution.

### Major changes

First stable release of advanced implementation of plant


