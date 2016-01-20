---
title: "plant: A package for modelling forest trait ecology & evolution: _Introduction_"
---

This document provides a number of worked examples that replicate
the figures in the paper.  The code is not identical to the paper
version because through this document we try for a bit more
exposition at the cost of more global variables, etc.

# Design

The `plant` package is designed around a set of nested components.

* At the highest level there is a metapopulation, represented by the
`SCM` object. (`SCM` denotes "Solver Characteristic Method".)

* The `SCM` object contains a single `Patch` object, but we use a
technique described in the "demography" document to demonstrate how
we can treat this patch as an infinite number of patches.

* The `Patch` object contains one or more `Species` objects.

* Each `Species` contains a number of of `Cohort` objects, each of
which contains a single `Plant` object representing the first
plant born in that cohort.

* Each `Plant` object has a `Strategy` which contains all the
physiological rules around growing a plant.  This `Strategy` is
shared across all members of a `Species`.

The document here works through these components in the reverse
order, starting with `Plant` (and really `Strategy`) objects, and
finishing by showing how fitness emerges from the `SCM`.  The
physiological model is described in detail in a separate "physiology"
document, while the details about demographic calculations are described
in the "demography" document.
