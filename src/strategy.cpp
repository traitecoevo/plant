#include <tree2/strategy.h>
#include <RcppCommon.h> // NA_REAL

namespace tree2 {

// TODO: There is some fairly major re-plumbing needed here; we need
// to separate out hyperparameters from the ones the model cares
// about, and possibly define a post-parameter setting hook so that
// other intermediates can be generated.

// TODO: Document consistent argument order: l, b, s, h, r
// TODO: Document ordering of different types of variables (size
// before physiology, before compound things?)

Strategy::Strategy() {
  // * Core traits - default values
  lma_0  = 0.1978791;  // leaf mas per area
  rho_0  = 608;        // wood density
  hmat_0 = 16.5958691; // Height at maturation
  s_0    = 3.8e-5;  // Seed size
  n_area_0 = 1.87e-3; // Leaf nitrogen per area (= Plant::v) [kg / m2]

  // * To start with set actiual values to default
  // TODO: This needs replacing with hyperparameters?
  lma    = lma_0;
  rho    = rho_0;
  hmat   = hmat_0;
  s      = s_0;
  n_area = n_area_0;

  // * Individual allometry
  // Canopy shape parameter (extra calculation here later)
  eta = 12;
  // ratio leaf area to sapwood area
  theta  = 4669;
  // Height - leaf mass scaling
  a1     = 5.44;
  B1     = 0.306;
  // Root - leaf scaling
  a3     = 0.07;
  // Ratio of bark area : sapwood area
  b      = 0.17;

  // * Production
  // Ratio of leaf dark respiration to leaf nitrogen mass
  // [mol CO2 / kgN / yr] (6.66e-4 * (365*24*60*60))
  c_Rl   = 2.1e4;
  // Root respiration per mass [mol CO2 / kg / yr]
  c_Rr   = 217;
  // Sapwood respiration per stem volume [mol CO2 / m3 / yr]
  c_Rs   = 4012;
  // Bark respiration per stem volume [mol CO2 / m3 / yr]
  // (note, new since paper -- see respiration calculation)
  c_Rb   = 2 * c_Rs;
  // Carbon conversion parameter
  Y      = 0.7;
  // Constant converting assimilated CO2 to dry mass [kg / mol]
  // (12E-3 / 0.49)
  c_bio  = 2.45e-2;
  // Leaf turnover - LMA scaling
  k_l0   =  0.4565855;
  B4     = 1.71;
  // Bark turnover
  k_b    = 0.2;
  // Sapwood turnover
  k_s0    = 0.2;
  B5      = 0;
  // Root turnover
  k_r    = 1.0;
  // Parameters of the hyperbola for annual LRC
  c_p1   = 150.36;
  c_p2   = 0.19;

  // * Seed production
  // Accessory cost of reproduction, kg for seed
  // with mass s_0
  // TODO: Deal with hyperparameters here?
  c_acc  = 3.0 * s_0;
  // Scaling of seed accessory costs with seed mass
  B7     = 1.0;

  // Maximum alloction to reproduction
  c_r1   = 1.0;
  // Size range across which individuals mature
  c_r2   = 50;

  // * Mortality parameters
  // Parameter for seedling survival
  c_s0    = 0.1;
  // Baseline for intrinsic mortality
  c_d0    = 0.01;
  // Coefficient for wood density in mortality function
  c_d1    = 0.0;
  // Coefficient for height density in mortality function
  B6    = 0.0;
  // Baseline rate for growth-related mortality
  c_d2    = 5.5;
  // Risk coefficient for dry mass production (per area)
  c_d3    = 20.0;

  // Will get computed properly by Plant.
  height_0    = NA_REAL;

  // NOTE: These need setting at some point
  // TODO: These possibly need setting as hyperparameters?
  eta_c = NA_REAL;
  k_l = NA_REAL;
  k_s = NA_REAL;
}

// [eqn 4] Sapwood area
double Strategy::sapwood_area(double leaf_area) const {
  return leaf_area / theta;
}

// [eqn 13] Total maintenance respiration
//
// (NOTE that there is a reparametrisation here relative to the paper
// -- c_Rb is defined (new) as 2*c_Rs, wheras the paper assumes a
// fixed multiplication by 2)
//
// NOTE: In contrast with EBT, we do not normalise by Y*c_bio.
//
// TODO: split into components?
double Strategy::respiration(double leaf_area, double sapwood_mass,
                             double bark_mass, double root_mass) const {
  return
    c_Rl * leaf_area * n_area +
    c_Rs * sapwood_mass / rho +
    c_Rb * bark_mass    / rho +
    c_Rr * root_mass;
}

// [eqn 14] Total turnover
//
// (NOTE: `k_l` is (a_4*\phi)^{b_4} in [eqn 14], and is computed by
// `prepare_strategy`).
// TODO[SPLIT]: done for sapwood but not for others
double Strategy::turnover(double leaf_mass, double bark_mass,
                          double sapwood_mass, double root_mass) const {
  return
    leaf_mass * k_l  +
    bark_mass * k_b  +
    sapwood_turnover(sapwood_mass) +
    root_mass * k_r;
}

// [eqn 15] Net production
//
// NOTE: Translation of variable names from the EBT.  everything
// before the minus sign is EBT's N, our `net_production` is EBT's P.
double Strategy::net_production(double assimilation, double respiration,
                                double turnover) {
  return c_bio * Y * (assimilation - respiration) - turnover;
}

// [eqn 16] Fraction of production allocated to reproduction
double Strategy::reproduction_fraction(double height) const {
  return c_r1 / (1.0 + exp(c_r2 * (1.0 - height / hmat)));
}

// [eqn 17] Rate of offspring production
// TODO[HYPERPARAMETER]: s/s0; ideally s + constant
double Strategy::dfecundity_dt(double net_production,
                               double reproduction_fraction) const {
  return net_production * reproduction_fraction /
    (s + c_acc * pow(s / s_0, B7));
}

// [eqn 18] Fraction of mass growth that is leaves (see doc/details.md
// for derivation).
double Strategy::leaf_mass_fraction(double leaf_area) const {
  return 1.0/(1.0
              + dsapwood_mass_dleaf_mass(leaf_area)
              + dbark_mass_dleaf_mass(leaf_area)
              + droot_mass_dleaf_mass(leaf_area));
}
double Strategy::leaf_area_fraction(double leaf_area) const {
  return leaf_mass_fraction(leaf_area) / lma;
}


// TODO: Ordering below here needs working on, probably as @dfalster
// does equation documentation?
double Strategy::dheight_dleaf_area(double leaf_area) const {
  return a1 * B1 * pow(leaf_area, B1 - 1);
}

// Mass of stem needed for new unit mass leaf, d m_s / d m_l
double Strategy::dsapwood_mass_dleaf_mass(double leaf_area) const {
  return rho * eta_c * a1 / (theta * lma) * (B1 + 1.0) * pow(leaf_area, B1);
}

// Mass of bark needed for new unit mass leaf, d m_b / d m_l
double Strategy::dbark_mass_dleaf_mass(double leaf_area) const {
  return b * dsapwood_mass_dleaf_mass(leaf_area);
}

// Mass of root needed for new unit mass leaf, d m_r / d m_l
double Strategy::droot_mass_dleaf_mass(double /* leaf_area */) const {
  return a3 / lma;
}

// Growth rate of basal diameter per unit time
double Strategy::dbasal_diam_dbasal_area(double leaf_area,
                                         double heartwood_area) const {
  return pow(M_PI / basal_area(leaf_area, heartwood_area), 0.5);
}

// Growth rate of leaf area per unit time
double Strategy::dleaf_area_dt(double leaf_mass_growth_rate) const {
  return leaf_mass_growth_rate / lma;
}

// Growth rate of spawood area at base per unit time
double Strategy::dsapwood_area_dt(double leaf_mass_growth_rate) const {
  return dleaf_area_dt(leaf_mass_growth_rate) / theta;
}

// TODO: @dfalster - should this not take leaf_mass_growth_rate like
// the others?
double Strategy::dheartwood_area_dt(double leaf_area) const {
  return k_s0 * sapwood_area(leaf_area);
}

// Growth rate of bark area at base per unit time
// TODO: this seems possible inefficient, probably does not matter
double Strategy::dbark_area_dt(double leaf_mass_growth_rate) const {
  return b * dsapwood_area_dt(leaf_mass_growth_rate);
}

// Growth rate of stem basal area per unit time
double Strategy::dbasal_area_dt(double leaf_area,
                                double leaf_mass_growth_rate) const {
  return dsapwood_area_dt(leaf_mass_growth_rate) +
    dbark_area_dt(leaf_mass_growth_rate) +
    dheartwood_area_dt(leaf_area);
}

// Growth rate of basal diameter per unit time
double Strategy::dbasal_diam_dt(double leaf_area,
                                double heartwood_area,
                                double leaf_mass_growth_rate) const {
  return dbasal_diam_dbasal_area(leaf_area, heartwood_area) *
    dbasal_area_dt(leaf_area, leaf_mass_growth_rate);
}

// TODO: Passing in leaf *area* but d (leaf *mass*) / dt, which does
// not seem ideal.
double Strategy::droot_mass_dt(double leaf_area,
                               double leaf_mass_growth_rate) const {
  return leaf_mass_growth_rate * droot_mass_dleaf_mass(leaf_area);
}

double Strategy::dlive_mass_dt(double reproduction_fraction,
                               double net_production) const {
  return (1 - reproduction_fraction) * net_production;
}

// TODO: Change top two to use dlive_mass_dt
double Strategy::dtotal_mass_dt(double reproduction_fraction,
                                double net_production,
                                double dheartwood_mass_dt) const {
  return dlive_mass_dt(reproduction_fraction, net_production) +
    dheartwood_mass_dt;
}

// TODO: Change top two to use dlive_mass_dt
// TODO: Do we not track root mass change?
double Strategy::dabove_ground_mass_dt(double leaf_area,
                                       double reproduction_fraction,
                                       double net_production,
                                       double dheartwood_mass_dt,
                                       double leaf_mass_growth_rate) const {
  const double droot_mass_dt =
    leaf_mass_growth_rate * droot_mass_dleaf_mass(leaf_area);
  return dtotal_mass_dt(reproduction_fraction, net_production,
                        dheartwood_mass_dt) - droot_mass_dt;
}

double Strategy::bark_area(double leaf_area) const {
  return b * sapwood_area(leaf_area);
}

double Strategy::basal_area(double leaf_area, double heartwood_area) const {
  return heartwood_area + bark_area(leaf_area) + sapwood_area(leaf_area);
}

double Strategy::dheartwood_mass_dt(double sapwood_mass) const {
  return sapwood_turnover(sapwood_mass);
}

double Strategy::sapwood_turnover(double sapwood_mass) const {
  return sapwood_mass * k_s;
}

double Strategy::height_given_leaf_mass(double leaf_mass) const {
  return a1 * pow(leaf_mass / lma, B1);
}

}
