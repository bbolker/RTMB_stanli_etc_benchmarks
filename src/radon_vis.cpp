// posteriordb radon_variable_intercept_slope_noncentered: negative log posterior
// on Stan's unconstrained scale (the sigmas are log-transformed)
#include <TMB.hpp>

template <class Type>
Type objective_function<Type>::operator()() {
  DATA_IVECTOR(county_idx);  // 1-based, as in the Stan data
  DATA_VECTOR(floor_measure);
  DATA_VECTOR(log_radon);
  PARAMETER(sigma_y);
  PARAMETER(sigma_alpha);
  PARAMETER(sigma_beta);
  PARAMETER_VECTOR(alpha_raw);
  PARAMETER_VECTOR(beta_raw);
  PARAMETER(mu_alpha);
  PARAMETER(mu_beta);

  Type sy = exp(sigma_y), sa = exp(sigma_alpha), sb = exp(sigma_beta);
  vector<Type> alpha = mu_alpha + sa * alpha_raw;
  vector<Type> beta = mu_beta + sb * beta_raw;

  Type lp = dnorm(sy, Type(0), Type(1), true) + dnorm(sb, Type(0), Type(1), true) +
            dnorm(sa, Type(0), Type(1), true) + dnorm(mu_alpha, Type(0), Type(10), true) +
            dnorm(mu_beta, Type(0), Type(10), true) +
            dnorm(alpha_raw, Type(0), Type(1), true).sum() +
            dnorm(beta_raw, Type(0), Type(1), true).sum();
  for (int i = 0; i < log_radon.size(); i++) {
    int j = county_idx(i) - 1;
    lp += dnorm(log_radon(i), alpha(j) + floor_measure(i) * beta(j), sy, true);
  }
  lp += sigma_y + sigma_alpha + sigma_beta;  // log-Jacobians of exp()
  return -lp;
}
