// -*-c++-*-
#ifndef PLANT_GRADIENT_SCM_HARVEST_H_
#define PLANT_GRADIENT_SCM_HARVEST_H_

// Strategy-agnostic helpers shared by the emergent-gradient engines (#472 scope B,
// refactor+optimize phase). These collapse boilerplate that was hand-copied across
// ff16_emergent.cpp / tf24_emergent.cpp / tf24f_emergent.cpp once the harvest moved
// into C++ (the native-SCM read + birth-step computation), plus the census
// height-trapezium reduction. Templated so the same code serves every strategy /
// scalar type; header-only so each engine TU instantiates what it needs.

#include <vector>
#include <cstddef>
#include <cmath>

namespace plant {
namespace gradient {

// Recover the constant birth rate from a run if not supplied (given < 0): the SCM's
// offspring_production / net_reproduction_ratio for the species (both per-species).
// `scm` is anything with operator-> to the SCM (a raw SCM* or an RcppR6<SCM> wrapper).
template <typename Scm>
double recover_birth_rate(Scm& scm, std::size_t species, double given) {
  if (given >= 0.0) return given;
  return scm->offspring_production()[species] /
         scm->net_reproduction_ratios()[species];
}

// Per-cohort birth steps: the (0-based) step time nearest each node introduction time.
// Mirrors the which.min(abs(sh - node_time)) the R harvests used.
template <typename Patch>
std::vector<int> birth_steps(const Patch& patch, std::size_t species) {
  const std::vector<double>& sh = patch.step_history;
  const std::vector<double> nt = patch.at_species(species).node_times();
  std::vector<int> birth(nt.size());
  for (std::size_t i = 0; i < nt.size(); ++i) {
    std::size_t best = 0;
    double bd = std::abs(sh[0] - nt[i]);
    for (std::size_t k = 1; k < sh.size(); ++k) {
      const double d = std::abs(sh[k] - nt[i]);
      if (d < bd) { bd = d; best = k; }
    }
    birth[i] = static_cast<int>(best);
  }
  return birth;
}

// The offspring-tape harvest computed natively from a live Patch: the per-cohort
// birth steps, node-spacing trapezoid weights (so offspring_production = sum_i tw_i *
// offspring_i), per-RK-stage patch survival, and survival-at-birth. Strategy-agnostic
// (every Patch exposes step_history / at_species / r_pr_survival / environment_history);
// shared by FF16's build_frozen_scm and the TF24/TF24f native offspring entries so the
// O(stand) survival loop runs once in C++ instead of being rebuilt in R. ppsurv is
// returned row-major [N*6] (caller wraps into the engine's matrix form).
struct OffspringWeights {
  std::size_t N;                // = environment_history.size() (Cash-Karp steps)
  std::vector<int> birth;       // per-cohort birth step (0-based)
  std::vector<double> tw;       // trapezoid node weight * patch_density * S_D * birth_rate
  std::vector<double> ppsab;    // pr_patch_survival_at_birth (per cohort)
  std::vector<double> ppsurv;   // [N*6] row-major: pr_survival at sh[k] + ah[s]*h_k
};
template <typename Patch>
OffspringWeights offspring_weights(const Patch& patch, std::size_t species,
                                   double S_D, double birth_rate) {
  OffspringWeights w;
  const std::vector<double>& sh = patch.step_history;
  const auto& sp = patch.at_species(species);
  const std::vector<double> nt    = sp.node_times();
  const std::vector<double> pdens = sp.r_patch_densities();
  w.ppsab = sp.r_pr_patch_survival_at_birth();
  const std::size_t nC = nt.size();
  w.birth = birth_steps(patch, species);

  std::vector<double> tcoef(nC, 0.0);
  if (nC >= 2) {
    tcoef[0] = 0.5 * (nt[1] - nt[0]);
    tcoef[nC - 1] = 0.5 * (nt[nC - 1] - nt[nC - 2]);
    for (std::size_t i = 1; i + 1 < nC; ++i) tcoef[i] = 0.5 * (nt[i + 1] - nt[i - 1]);
  }
  w.tw.resize(nC);
  for (std::size_t i = 0; i < nC; ++i) w.tw[i] = tcoef[i] * pdens[i] * S_D * birth_rate;

  w.N = patch.environment_history.size();
  const double ah[6] = {0.0, 0.2, 0.3, 0.6, 1.0, 0.875};
  w.ppsurv.assign(w.N * 6, 0.0);
  for (std::size_t k = 0; k < w.N; ++k) {
    const double hN = sh[k + 1] - sh[k];
    for (int st = 0; st < 6; ++st)
      w.ppsurv[k * 6 + st] = patch.r_pr_survival(sh[k] + ah[st] * hN);
  }
  return w;
}

// Census height-trapezium over the replayed cohorts, descending height + the
// pending-seed tail term (mirrors Species::compute_competition / Patch census
// integral). Generic over the cohort container: `geth(i)` / `getdens(i)` / `getmhw(i)`
// extract the per-cohort height / number density / heartwood mass at index i, and
// `psi(h, dens, mhw)` is the metric kernel. `ord` indexes cohorts in descending
// height; (h0, dens_new) are the pending-seed height / density for the ground tail.
// Identical arithmetic to the hand-rolled reductions it replaces (bit-for-bit).
template <typename S, typename GetH, typename GetDens, typename GetMhw, typename Psi>
S census_trapezium(std::size_t nC, const std::vector<std::size_t>& ord, S h0, S dens_new,
                   GetH geth, GetDens getdens, GetMhw getmhw, Psi psi) {
  std::vector<S> phi(nC);
  for (std::size_t i = 0; i < nC; ++i) phi[i] = psi(geth(i), getdens(i), getmhw(i));
  S J = S(0.0);
  for (std::size_t j = 0; j + 1 < nC; ++j) {
    const std::size_t a = ord[j], b = ord[j + 1];
    J = J + S(0.5) * (geth(a) - geth(b)) * (phi[a] + phi[b]);
  }
  if (nC > 0) {
    const std::size_t last = ord[nC - 1];
    S phi_new = psi(h0, dens_new, S(0.0));
    J = J + S(0.5) * (geth(last) - h0) * (phi[last] + phi_new);
  }
  return J;
}

} // namespace gradient
} // namespace plant

#endif
