# R version of the growth model

This is my interface to Daniel's *R* version of the growth model, used
to cross check the C++ version implemented by this package.

Daniel's version makes extensive use of global variables for both
parameters and state.  I don't want to rewrite the model to change
this, partly because that introduces the possibility of changes that
alter the results, but partly because this code is only going to be
used for testing.

The workaround is that rather than keeping the files in the `R/`
directory, they are here in the `inst/` directory.  There is a
function in the package `make_reference_plant` that loads the contents of
these files into a new environment and makes an object that can be
manipulated and have parameters set.  This prevents any pollution of
the package environment.

The function `make_reference_plant` loads the code within the `R/`
directory.  The files `R/plots.R`, `R/main.R` and `R/main 2.R` are not
currently used.

The file `extra.R` (in *this* directory) is additional code added
during development of the C++ version of the model.

This code was emailed by Daniel (daniel.falster, mq.edu.au) to me
(rich.fitzjohn, gmail.com) on 2012/10/09 13:29 AEST (in
`Archive.zip`).
