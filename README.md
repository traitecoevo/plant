# The TRait Ecology and Evolution (TREE) model

[![Build Status](https://travis-ci.org/traitecoevo/tree.png?branch=master)](https://travis-ci.org/traitecoevo/tree)

## Installation

Installation requires a C++11 compatible C compiler (OSX >= 10.10/Yosemite satisfies this, plus we've had success on Ubuntu 12.4 and 14.4).

**Option 1, using [`drat`](https://github.com/eddelbuettel/drat)**

(install `drat` with `install.packages("drat")`)

```r
drat:::add("traitecoevo")
install.packages("tree")
```

(versions of R before 3.2.0 may require `install.packages("tree", type="source")` here)

**Option 2, using `devtools::install_github`**

```r
devtools::install_github("traitecoevo/tree")
```
