// -*-c++-*-
#ifndef PLANT_FP_CONTRACT_H_
#define PLANT_FP_CONTRACT_H_

// No fused multiply-add. A forward pass runs at double and a sweep records the
// same arithmetic through the active scalar's operators, and the two must agree
// to the last bit: a sweep differentiates the run it replays. Clang contracts
// a * b + c into one fma in double code by default, and never through an
// operator call, so the two can differ by an ulp -- which a leaf solve then
// amplifies to ~1e-8 (the state-Jacobian ladder, macOS only, after #700).
//
// ⚠️ THIS MUST PRECEDE EVERY PLANT AND ODELIA HEADER IN A TRANSLATION UNIT.
// Clang fixes a function's contraction where the function is DEFINED, so a
// template defined above this line keeps the default. plant.h includes it first,
// which covers the generated RcppExports.cpp and RcppR6.cpp; every other source
// file includes it as its first line. A pragma, not -ffp-contract=off in
// src/Makevars, because R CMD check reports any -f flag in PKG_CXXFLAGS as a
// non-portable WARNING.
//
// GCC has no equivalent for C++. It contracts only where the target has fma, so
// x86-64 built for the generic target (CRAN, the CI runners) is unaffected;
// aarch64 Linux under GCC is not covered.
#if defined(__clang__)
#pragma clang fp contract(off)
#endif

#endif
