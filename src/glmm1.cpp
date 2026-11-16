// posteriordb GLMM1_model: negative log posterior on Stan's unconstrained scale
#include <TMB.hpp>
#include "stan_transforms.hpp"

template <class Type>
Type objective_function<Type>::operator()() {
  DATA_VECTOR(obs);
  DATA_IVECTOR(obssite);  // 1-based, as in the Stan data
  PARAMETER_VECTOR(alpha);
  PARAMETER(mu_alpha);
  PARAMETER(sd_alpha);

  Type sd = lub_trans(sd_alpha, Type(0), Type(5));
  Type lp = dnorm(alpha, mu_alpha, sd, true).sum() + dnorm(mu_alpha, Type(0), Type(10), true);
  for (int i = 0; i < obs.size(); i++) lp += dpois(obs(i), exp(alpha(obssite(i) - 1)), true);
  lp += lub_logjac(sd_alpha, Type(0), Type(5));
  return -lp;
}
