// posteriordb GLMM_Poisson_model: negative log posterior on Stan's unconstrained scale
#include <TMB.hpp>
#include "stan_transforms.hpp"

template <class Type>
Type objective_function<Type>::operator()() {
  DATA_VECTOR(C);
  DATA_VECTOR(year);
  PARAMETER(alpha);
  PARAMETER(beta1);
  PARAMETER(beta2);
  PARAMETER(beta3);
  PARAMETER_VECTOR(eps);
  PARAMETER(sigma);

  Type a  = lub_trans(alpha, Type(-20), Type(20));
  Type b1 = lub_trans(beta1, Type(-10), Type(10));
  Type b2 = lub_trans(beta2, Type(-10), Type(10));
  Type b3 = lub_trans(beta3, Type(-10), Type(10));
  Type s  = lub_trans(sigma, Type(0), Type(5));
  vector<Type> eta = a + b1 * year + b2 * year * year + b3 * year * year * year + eps;

  Type lp = dpois(C, exp(eta), true).sum() + dnorm(eps, Type(0), s, true).sum();
  lp += lub_logjac(alpha, Type(-20), Type(20)) + lub_logjac(beta1, Type(-10), Type(10)) +
        lub_logjac(beta2, Type(-10), Type(10)) + lub_logjac(beta3, Type(-10), Type(10)) +
        lub_logjac(sigma, Type(0), Type(5));
  return -lp;
}
