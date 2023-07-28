#include <TMB.hpp>
#include "utils.h"

// using namespace density;
// using namespace Eigen;

//distribution enumeration
enum valid_family {
  gaussian_family = 0,
  poisson_family = 1,
};

enum valid_link {
  identity_link = 0,
  log_link = 1
};

// Inverse link function
template<class Type>
Type inv_link(Type eta, int link) {
  Type ans;
  switch (link) {
  case identity_link:
    ans = eta;
    break;
  case log_link:
    ans = exp(eta);
    break;
  default:
    error("Link not implemented!");
  } // End switch
  return ans;
}


//objective function -----------------------------------------------------------

template<class Type>
Type objective_function<Type>::operator() () {

  //data
  DATA_VECTOR(y);
  DATA_MATRIX(X);
  DATA_VECTOR(xi);
  DATA_STRUCT(term, snglmm::term_t);
  DATA_SPARSE_MATRIX(Z);
  DATA_MATRIX(D);
  DATA_MATRIX(W);

  DATA_INTEGER(link);
  DATA_INTEGER(family);

  //parameter
  PARAMETER_VECTOR(b);
  PARAMETER(log_sigma);
  PARAMETER_VECTOR(log_theta);
  PARAMETER_VECTOR(v);
  PARAMETER_VECTOR(u);
  PARAMETER(log_phi);
  PARAMETER(log_lambda);

  //transformed parameters
  Type sigma = exp(log_sigma);
  Type lambda = exp(log_lambda);
  Type phi = exp(log_phi);
  vector<Type> theta = exp(log_theta);

  //negative log likelihood
  Type nll, jnll = 0;

  //non-spatial random effects
  jnll += snglmm::mnll(v, log_theta, term);
  vector<Type> eta = Z * v;

  //spatial random effects
  if (D.sum() != 0) {
    jnll += snglmm::snll(u, log_phi, log_lambda, D, W);
    eta += u;
  }

  //fixed effects (xi is an offset term)
  eta += X * b + xi;
  vector<Type> mu(eta.size());

  for (int i = 0; i < eta.size(); i++)
    mu(i) = inv_link(eta(i), link);

  for (int i = 0; i < eta.size(); i++) {
    switch (family) {
    case gaussian_family:
      nll = dnorm(y(i), mu(i), sigma, true);
      break;
    case poisson_family:
      nll = dpois(y(i), mu(i), true);
      break;
    default:
      error("Family not implemented!");
    }//switch end
    jnll -= nll;
  }

  vector<matrix<Type> > cor(term.size());
  for(int i = 0; i < term.size(); i++){
    cor(i) = term(i).R;
  }

  //reporting ------------------------------------------------------------------

  REPORT(cor);
  REPORT(eta);
  REPORT(v);
  REPORT(u);

  ADREPORT(b);
  ADREPORT(sigma);
  ADREPORT(theta);
  ADREPORT(lambda);
  ADREPORT(phi);

  return jnll;
}
