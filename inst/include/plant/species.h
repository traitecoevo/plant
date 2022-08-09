// -*-c++-*-
#ifndef SPECIES
#define SPECIES

#include <vector>
#include <plant/util.h>
#include <plant/environment.h>
#include <plant/ode_interface.h>
#include <plant/node.h>
#include <plant/extrinsic_drivers.h>

namespace plant {

// This is purely for running the deterministic model.

template <typename T, typename E>
class Species {
public:
  typedef T         strategy_type;
  typedef E         environment_type;
  typedef Individual<T,E>  individual_type;
  typedef Node<T,E> node_type;
  typedef typename strategy_type::ptr strategy_type_ptr;
  Species(strategy_type s);

  size_t size() const;
  void clear();
  void introduce_new_node();

  double height_max() const;
  double compute_competition(double height) const;
  void compute_rates(const environment_type& environment, double pr_patch_survival, double birth_rate);
  std::vector<double> net_reproduction_ratio_by_node() const;

  // * ODE interface
  // NOTE: We are a time-independent model here so no need to pass
  // time in as an argument.  All the bits involving time are taken
  // care of by Environment for us.
  size_t ode_size() const;
  size_t aux_size() const;
  size_t strategy_aux_size() const;
  std::vector<std::string> aux_names() const;

  void resize_consumption_rates(int i);
  double consumption_rate(int i) const;
  std::vector<double> consumption_rate_by_node_rev(int i) const;

  ode::const_iterator set_ode_state(ode::const_iterator it);
  ode::iterator       ode_state(ode::iterator it) const;
  ode::iterator       ode_rates(ode::iterator it) const;
  ode::iterator       ode_aux(ode::iterator it) const;

  // * R interface
  std::vector<double> r_heights() const;
  std::vector<double> r_heights_rev() const;
  void r_set_heights(std::vector<double> heights);
  const node_type& r_new_node() const {return new_node;}
  std::vector<node_type> r_nodes() const {return nodes;}
  const node_type& r_node_at(util::index idx) const {
    return nodes[idx.check_bounds(size())];
  }

  // Do this with set_ode_state, using an iterator?
  /* double state(int i) const { return vars.state(i); } */

  /* double rate(int i) const { return vars.rate(i); } */

  /* void set_state(int i, double v) { */
  /*   vars.set_state(i, v); */
  /* } */


  // These are used to determine the degree of node refinement.
  std::vector<double> r_competition_effects() const;
  std::vector<double> r_competition_effects_error(double scal) const;

  // This is just kind of useful
  std::vector<double> r_log_densities() const;
  const ExtrinsicDrivers& extrinsic_drivers() const {return strategy->extrinsic_drivers;}

private:
  const Control& control() const {return strategy->get_control();}
  strategy_type_ptr strategy;
  node_type new_node;
  std::vector<node_type> nodes;

  typedef typename std::vector<node_type>::iterator nodes_iterator;
  typedef typename std::vector<node_type>::const_iterator nodes_const_iterator;
};

template <typename T, typename E>
Species<T,E>::Species(strategy_type s)
  : strategy(make_strategy_ptr(s)),
    new_node(strategy) {
}

template <typename T, typename E>
size_t Species<T,E>::size() const {
  return nodes.size();
}

template <typename T, typename E>
void Species<T,E>::clear() {
  nodes.clear();
  // Reset the new_node to a blank new_node, too.
  new_node = node_type(strategy);
}

template <typename T, typename E>
void Species<T,E>::introduce_new_node() {
  nodes.push_back(new_node);
  // TODO: Should the new_node be recomputed here?
}

// If a species contains no individuals, we return the height of a
// seed of the species.  Otherwise we return the height of the largest
// individual (always the first in the list) which will be at least
// tall as a seed.
template <typename T, typename E>
double Species<T,E>::height_max() const {
  return nodes.empty() ? new_node.height() : nodes.front().height();
}

// Because of nodes are always ordered from largest to smallest, we
// need not continue down the list once the leaf area above a certain
// height is zero, because it will be zero for all nodes further down
// the list.
//
// NOTE: This is simply performing numerical integration,  via the
// trapezium rule, of the compute_competition with respect to plant
// height.  You'd think that this would be nicer to do in terms of a
// call to an external trapezium integration function, but building
// and discarding the intermediate storage ends up being a nontrivial
// cost.  A more general iterator version might be possible, but with
// the fiddliness around the boundary conditions that won't likely be
// useful.
//
// NOTE: In the cases where there is no individuals, we return 0 for
// all heights.  The integral is not defined, but an empty light
// environment seems appropriate.
//
// NOTE: A similar early-exit condition to the Plant version is used;
// once the lower bound of the trazpeium is zero, we stop including
// individuals.  Working with the boundary node is tricky here,
// because we might need to include that, too: always in the case of a
// single node (needed to be the second half of the trapezium) and
// also needed if the last looked at plant was still contributing to
// the integral).
template <typename T, typename E>
double Species<T,E>::compute_competition(double height) const {
  if (size() == 0 || height_max() < height) {
    return 0.0;
  }
  double tot = 0.0;
  nodes_const_iterator it = nodes.begin();
  double h1 = it->height(), f_h1 = it->compute_competition(height);

  // Loop over nodes
  for (++it; it != nodes.end(); ++it) {
    
    const double h0 = it->height(), f_h0 = it->compute_competition(height);
    if (!util::is_finite(f_h0)) {
      util::stop("Detected non-finite contribution");
    }

    // Integration
    tot += (h1 - h0) * (f_h1 + f_h0);
    // Upper point moves for next time:
    h1   = h0;
    f_h1 = f_h0;
    if (h0 < height) {
      break;
    }
  }

  if (size() == 1 || f_h1 > 0) {
    const double h0 = new_node.height(), f_h0 = new_node.compute_competition(height);
    tot += (h1 - h0) * (f_h1 + f_h0);
  }

  return tot / 2;
}

// NOTE: We should probably prefer to rescale when this is called
// through the ode stepper.
template <typename T, typename E>
void Species<T,E>::compute_rates(const E& environment, double pr_patch_survival, double birth_rate) {
    int counter_species = 0; 
  std::cout << "\ntime\t"<< environment.time; 

  for (auto& c : nodes) {
  counter_species += 1;
  std::cout << "\tnode\t" << counter_species; 
    c.compute_rates(environment, pr_patch_survival);
  }
  new_node.compute_initial_conditions(environment, pr_patch_survival, birth_rate);
}

template <typename T, typename E>
std::vector<double> Species<T,E>::net_reproduction_ratio_by_node() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(c.fecundity());
  }
  return ret;
}

template <typename T, typename E>
void Species<T,E>::resize_consumption_rates(int r) {
  new_node.resize_consumption_rates(r);
}

template <typename T, typename E>
double Species<T,E>::consumption_rate(int i) const {
  // can't determine density for one node
  if(size() < 2) {
    return 0.0;
  } else {
    // node heights are in descending order - we need ascending for integration
    return util::trapezium(r_heights_rev(), consumption_rate_by_node_rev(i));
  }
}

template <typename T, typename E>
std::vector<double> Species<T,E>::consumption_rate_by_node_rev(int i) const {
  std::vector<double> ret;
  ret.reserve(size());
  for(auto it = nodes.rbegin(); it != nodes.rend(); ++it) {
    ret.push_back(it->consumption_rate(i));
  }
  return ret;
}

template <typename T, typename E>
size_t Species<T,E>::ode_size() const {
  return size() * node_type::ode_size();
}

// bit clunky...
template <typename T, typename E>
size_t Species<T,E>::aux_size() const {
  return size() * strategy_aux_size();
}

// these 2 only really used in get_aux.h
template <typename T, typename E>
size_t Species<T,E>::strategy_aux_size() const {
  return strategy->aux_size();
}

template <typename T, typename E>
std::vector<std::string> Species<T,E>::aux_names() const {
  return strategy->aux_names();
}

template <typename T, typename E>
ode::const_iterator Species<T,E>::set_ode_state(ode::const_iterator it) {
  return ode::set_ode_state(nodes.begin(), nodes.end(), it);
}

template <typename T, typename E>
ode::iterator Species<T,E>::ode_state(ode::iterator it) const {
  return ode::ode_state(nodes.begin(), nodes.end(), it);
}

template <typename T, typename E>
ode::iterator Species<T,E>::ode_rates(ode::iterator it) const {
  return ode::ode_rates(nodes.begin(), nodes.end(), it);
}
//double sum_aux(int index) {}

template <typename T, typename E>
ode::iterator Species<T,E>::ode_aux(ode::iterator it) const {
  return ode::ode_aux(nodes.begin(), nodes.end(), it);
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(it->height());
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_heights_rev() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(it->height());
  }
  std::reverse(ret.begin(), ret.end());
  return ret;
}

template <typename T, typename E>
void Species<T,E>::r_set_heights(std::vector<double> heights) {
  util::check_length(heights.size(), size());
  if (!util::is_decreasing(heights.begin(), heights.end())) {
    util::stop("height must be decreasing (ties allowed)");
  }
  size_t i = 0;
  for (nodes_iterator it = nodes.begin(); it != nodes.end(); ++it, ++i) {
    it->individual.set_state("height", heights[i]);
  }
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_competition_effects() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (auto& c : nodes) {
    ret.push_back(c.competition_effect());
  }
  return ret;
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_competition_effects_error(double scal) const {
  return util::local_error_integration(r_heights(), r_competition_effects(), scal);
}

template <typename T, typename E>
std::vector<double> Species<T,E>::r_log_densities() const {
  std::vector<double> ret;
  ret.reserve(size());
  for (nodes_const_iterator it = nodes.begin();
       it != nodes.end(); ++it) {
    ret.push_back(it->get_log_density());
  }
  return ret;
}

}

#endif /* SPECIES */
