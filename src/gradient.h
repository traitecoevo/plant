// -*-c++-*-
#ifndef TREE_GRADIENT_H_
#define TREE_GRADIENT_H_

#include "functor.h"

namespace util {

double gradient_fd_forward(DFunctor *f, double x, double dx, double fx);
double gradient_fd_forward(DFunctor *f, double x, double dx);
double gradient_fd_centre(DFunctor *f, double x, double dx);
double gradient_richardson(DFunctor *f, double x);

namespace test {
double test_gradient(double x, double dx, bool centre,
		     std::vector<double> pars);
double test_gradient_richardson(double x, std::vector<double> pars);
}

}

#endif
