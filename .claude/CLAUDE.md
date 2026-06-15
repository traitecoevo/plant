# CLAUDE.md

This project's architecture and contributor guidance for AI agents lives in
[agents.md](../agents.md). **Read [agents.md](../agents.md) before making changes** —
especially anything touching the C++/R boundary, the RcppR6 interface, or
adding a new model.

Quick pointers:

- Two-layer package: C++ core (`inst/include/plant/`, `src/`) bridged to R
  (`R/`) via RcppR6 + Rcpp. Interface contract: `inst/RcppR6_classes.yml`.
- **Never hand-edit generated files** (`R/RcppR6.R`, `R/RcppExports.R`,
  `src/RcppR6.cpp`, `src/RcppExports.cpp`, `NAMESPACE`, `man/`).
- After interface changes run `make rebuild`; after C++-only changes
  `make compile`; for R-only changes use `devtools::load_all()`.
- Run tests with `make test` or `devtools::test()`.

See [agents.md](../agents.md) for the full architecture, build workflow, and the
process for adding a new strategy/environment model.
