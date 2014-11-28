#include <tree/functor.h>

#include <tree/util.h>

namespace util {

DFunctor::~DFunctor() {
}

double helper_functor(double x, void *data) {
  DFunctor *f = static_cast<DFunctor *>(data);
  return (*f)(x);
}

RFunctionWrapper::RFunctionWrapper(Rcpp::Function fun_)
  : fun(fun_) {
}

double RFunctionWrapper::target(double x) {
  return Rcpp::as<double>(fun(Rcpp::wrap(x)));
}

double RFunctionWrapper::operator()(double x) {
  return target(x);
}

namespace test {

std::vector<double> test_functor(std::vector<double> x, 
				 std::vector<double> pars) {
  util::check_length(pars.size(), 3);
  Quadratic obj(pars[0], pars[1], pars[2]);
  Functor<Quadratic, &Quadratic::mytarget> fun(&obj);
  
  std::vector<double> ret;
  for (std::vector<double>::iterator it = x.begin(); it != x.end(); ++it)
    ret.push_back(fun(*it));
  
  return ret;
}

}

}
