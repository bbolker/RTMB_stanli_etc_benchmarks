// posteriordb GLM_Binomial_model: negative log posterior (no constrained parameters)
#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()() {
  DATA_VECTOR(C);
  DATA_VECTOR(N);
  DATA_VECTOR(year);
  PARAMETER(alpha);
  PARAMETER(beta1);
  PARAMETER(beta2);

  vector<Type> eta = alpha + beta1 * year + beta2 * year * year;
  Type lp = dnorm(alpha, Type(0), Type(100), true) + dnorm(beta1, Type(0), Type(100), true) +
            dnorm(beta2, Type(0), Type(100), true);
  for (int i = 0; i < C.size(); i++) lp += dbinom_robust(C(i), N(i), eta(i), true);
  return -lp;
}
