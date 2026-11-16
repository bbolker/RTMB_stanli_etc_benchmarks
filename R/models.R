## Model definitions for the benchmarks: posteriordb Stan code and data, plus
## equivalent RTMB, TMB, glmmTMB and brms versions.
##
## The RTMB and TMB objectives are the *negative log posterior on Stan's
## unconstrained scale*, Jacobian included, with parameters in Stan's
## declaration order. So obj$fn()/obj$gr() can be compared point by point with
## cmdstanr's log_prob()/grad_log_prob() and stanli's log_prob_grad() (up to an
## additive constant: Stan's `~` drops normalizing constants), and
## tmbstan/StanEstimators sample exactly the same target as Stan. Parameter
## names follow the Stan code, but the values are unconstrained (e.g.
## `sigma_y` holds log(sigma_y)).
##
## glmmTMB and brms are the "typical user" versions: their own
## parameterizations and priors, so they are not part of the
## log-density check.
##
## To add a model, add an entry to `pdb_names` and to `model_specs`, and a
## matching src/<short>.cpp.

library(RTMB)

model_dir <- here::here("models")
src_dir <- here::here("src")

## short name -> posteriordb posterior name
pdb_names <- c(
  glm_pois  = "GLM_Poisson_Data-GLM_Poisson_model",
  glm_binom = "GLM_Binomial_data-GLM_Binomial_model",
  glmm_pois = "GLMM_Poisson_data-GLMM_Poisson_model",
  glmm1     = "GLMM_data-GLMM1_model",
  radon_vis = "radon_mn-radon_variable_intercept_slope_noncentered"
)

## Stan code and data, cached in models/<short>.{stan,json} so posteriordb is
## only contacted once.
get_pdb_model <- function(short) {
  stan_file <- file.path(model_dir, paste0(short, ".stan"))
  json_file <- file.path(model_dir, paste0(short, ".json"))
  if (!file.exists(stan_file) || !file.exists(json_file)) {
    dir.create(model_dir, showWarnings = FALSE)
    po <- posteriordb::posterior(pdb_names[[short]], posteriordb::pdb_github())
    writeLines(posteriordb::stan_code(po), stan_file)
    ## auto_unbox turns length-1 vectors into scalars, so a length-1 data
    ## *array* would be written as a scalar (none of the current models has one)
    jsonlite::write_json(posteriordb::get_data(po), json_file,
                         auto_unbox = TRUE, digits = NA)
  }
  list(stan_file = stan_file, json_file = json_file,
       data = jsonlite::read_json(json_file, simplifyVector = TRUE))
}

## Stan's transforms and log-Jacobians (work on numerics and RTMB advectors)
lub_trans  <- function(u, a, b) a + (b - a) / (1 + exp(-u))
lub_logjac <- function(u, a, b) log(b - a) - logspace_add(0 * u, -u) - logspace_add(0 * u, u)
## lower = 0: x = exp(u), log-Jacobian = u

## glmmTMB `priors` data frame. glmmTMB's normal prior on an SD parameter acts
## on log-SD; its gamma prior acts on the SD but adds no Jacobian, so when
## tmbstan samples log-SD, gamma(mean, shape = 2) is an exponential prior on
## the SD with scale mean/2. `fixef` skips the intercept, hence `intercept`.
## `disp` is a prior on the (single) dispersion parameter, e.g. the residual
## SD; it is given by index because glmmTMB (1.1.15.9000) can't match
## fixef_disp elements by name (and a blank coef there applies no prior).
glmmTMB_prior_df <- function(intercept = NULL, fixef = NULL, ranef_sd = NULL, disp = NULL) {
  rows <- list(
    if (!is.null(intercept)) c(intercept, "fixef", "(Intercept)"),
    if (!is.null(fixef)) c(fixef, "fixef", ""),
    if (!is.null(ranef_sd)) c(ranef_sd, "ranef_sd", ""),
    if (!is.null(disp)) c(disp, "fixef_disp", "1")
  )
  m <- do.call(rbind, rows)
  data.frame(prior = m[, 1], class = m[, 2], coef = m[, 3])
}
## exponential(scale 100) on the SD: almost flat, standing in for uniform(0, 5)
sd_flat <- "gamma(200,2)"
## exponential with the mean of a half-normal(0, 1) (0.8), standing in for it
sd_halfnorm1 <- "gamma(1.6,2)"

## Per-model pieces. `par(d)`: unconstrained starting values in Stan's
## declaration order; `lp(p, d)`: log posterior (unconstrained scale);
## `glmmTMB(d)`/`brms(d)`: arguments for those packages.
model_specs <- list(

  glm_pois = list(
    par = function(d) list(alpha = 0, beta1 = 0, beta2 = 0, beta3 = 0),
    lp = function(p, d) {
      a  <- lub_trans(p$alpha, -20, 20)
      b1 <- lub_trans(p$beta1, -10, 10)
      b2 <- lub_trans(p$beta2, -10, 10)
      b3 <- lub_trans(p$beta3, -10, 10)
      eta <- a + b1 * d$year + b2 * d$year^2 + b3 * d$year^3
      sum(dpois(d$C, exp(eta), log = TRUE)) +
        lub_logjac(p$alpha, -20, 20) + lub_logjac(p$beta1, -10, 10) +
        lub_logjac(p$beta2, -10, 10) + lub_logjac(p$beta3, -10, 10)
    },
    df = function(d) data.frame(C = d$C, year = d$year),
    glmmTMB = list(formula = C ~ year + I(year^2) + I(year^3), family = poisson),
    ## flat priors approximate the wide uniform priors
    glmmTMB_priors = NULL,
    brms = list(formula = C ~ year + I(year^2) + I(year^3), family = poisson,
                prior = NULL)
  ),

  glm_binom = list(
    par = function(d) list(alpha = 0, beta1 = 0, beta2 = 0),
    lp = function(p, d) {
      eta <- p$alpha + p$beta1 * d$year + p$beta2 * d$year^2
      dnorm(p$alpha, 0, 100, log = TRUE) + dnorm(p$beta1, 0, 100, log = TRUE) +
        dnorm(p$beta2, 0, 100, log = TRUE) +
        sum(dbinom_robust(d$C, d$N, eta, log = TRUE))
    },
    df = function(d) data.frame(C = d$C, N = d$N, year = d$year),
    glmmTMB = list(formula = cbind(C, N - C) ~ year + I(year^2), family = binomial),
    glmmTMB_priors = glmmTMB_prior_df(fixef = "normal(0,100)", intercept = "normal(0,100)"),
    brms = list(formula = C | trials(N) ~ year + I(year^2), family = "binomial",
                prior = c(brms::prior(normal(0, 100), class = b),
                          brms::prior(normal(0, 100), class = Intercept)))
  ),

  glmm_pois = list(
    par = function(d) list(alpha = 0, beta1 = 0, beta2 = 0, beta3 = 0,
                           eps = rep(0, d$n), sigma = 0),
    lp = function(p, d) {
      a  <- lub_trans(p$alpha, -20, 20)
      b1 <- lub_trans(p$beta1, -10, 10)
      b2 <- lub_trans(p$beta2, -10, 10)
      b3 <- lub_trans(p$beta3, -10, 10)
      sigma <- lub_trans(p$sigma, 0, 5)
      eta <- a + b1 * d$year + b2 * d$year^2 + b3 * d$year^3 + p$eps
      sum(dpois(d$C, exp(eta), log = TRUE)) +
        sum(dnorm(p$eps, 0, sigma, log = TRUE)) +
        lub_logjac(p$alpha, -20, 20) + lub_logjac(p$beta1, -10, 10) +
        lub_logjac(p$beta2, -10, 10) + lub_logjac(p$beta3, -10, 10) +
        lub_logjac(p$sigma, 0, 5)
    },
    df = function(d) data.frame(C = d$C, year = d$year, obs = factor(seq_len(d$n))),
    glmmTMB = list(formula = C ~ year + I(year^2) + I(year^3) + (1 | obs),
                   family = poisson),
    glmmTMB_priors = glmmTMB_prior_df(ranef_sd = sd_flat),
    brms = list(formula = C ~ year + I(year^2) + I(year^3) + (1 | obs),
                family = poisson, prior = NULL)
  ),

  glmm1 = list(
    ## no year effect in the reference model (log_lambda doesn't vary by year)
    par = function(d) list(alpha = rep(0, d$nsite), mu_alpha = 0, sd_alpha = 0),
    lp = function(p, d) {
      sd_alpha <- lub_trans(p$sd_alpha, 0, 5)
      sum(dnorm(p$alpha, p$mu_alpha, sd_alpha, log = TRUE)) +
        dnorm(p$mu_alpha, 0, 10, log = TRUE) +
        sum(dpois(d$obs, exp(p$alpha[d$obssite]), log = TRUE)) +
        lub_logjac(p$sd_alpha, 0, 5)
    },
    df = function(d) data.frame(obs = d$obs, site = factor(d$obssite)),
    glmmTMB = list(formula = obs ~ 1 + (1 | site), family = poisson),
    glmmTMB_priors = glmmTMB_prior_df(intercept = "normal(0,10)", ranef_sd = sd_flat),
    brms = list(formula = obs ~ 1 + (1 | site), family = poisson,
                prior = brms::prior(normal(0, 10), class = Intercept))
  ),

  radon_vis = list(
    par = function(d) list(sigma_y = 0, sigma_alpha = 0, sigma_beta = 0,
                           alpha_raw = rep(0, d$J), beta_raw = rep(0, d$J),
                           mu_alpha = 0, mu_beta = 0),
    lp = function(p, d) {
      sigma_y <- exp(p$sigma_y)
      sigma_alpha <- exp(p$sigma_alpha)
      sigma_beta <- exp(p$sigma_beta)
      alpha <- p$mu_alpha + sigma_alpha * p$alpha_raw
      beta <- p$mu_beta + sigma_beta * p$beta_raw
      cid <- d$county_idx
      mu <- alpha[cid] + d$floor_measure * beta[cid]
      dnorm(sigma_y, 0, 1, log = TRUE) + dnorm(sigma_beta, 0, 1, log = TRUE) +
        dnorm(sigma_alpha, 0, 1, log = TRUE) +
        dnorm(p$mu_alpha, 0, 10, log = TRUE) + dnorm(p$mu_beta, 0, 10, log = TRUE) +
        sum(dnorm(p$alpha_raw, 0, 1, log = TRUE)) +
        sum(dnorm(p$beta_raw, 0, 1, log = TRUE)) +
        sum(dnorm(d$log_radon, mu, sigma_y, log = TRUE)) +
        p$sigma_y + p$sigma_alpha + p$sigma_beta
    },
    df = function(d) data.frame(log_radon = d$log_radon, floor = d$floor_measure,
                                county = factor(d$county_idx)),
    glmmTMB = list(formula = log_radon ~ floor + diag(1 + floor | county),
                   family = gaussian),
    glmmTMB_priors = glmmTMB_prior_df(intercept = "normal(0,10)", fixef = "normal(0,10)",
                                      ranef_sd = sd_halfnorm1, disp = sd_halfnorm1),
    ## brms has no diag(); `||` gives the same uncorrelated intercept and slope
    brms = list(formula = log_radon ~ floor + (1 + floor || county),
                family = gaussian,
                prior = c(brms::prior(normal(0, 10), class = Intercept),
                          brms::prior(normal(0, 10), class = b),
                          brms::prior(normal(0, 1), class = sd),
                          brms::prior(normal(0, 1), class = sigma)))
  )
)

model_names <- names(pdb_names)
stopifnot(identical(sort(names(model_specs)), sort(model_names)))

## ---- builders ----

## RTMB: obj$fn is the negative log posterior (unconstrained scale)
make_rtmb <- function(short, d = get_pdb_model(short)$data) {
  spec <- model_specs[[short]]
  MakeADFun(function(p) -spec$lp(p, d), spec$par(d), silent = TRUE)
}

## TMB: compile src/<short>.cpp (if needed) and build the same objective.
## TMB's DATA_* macros want doubles, so integer data are converted.
compile_tmb <- function(short, ...) {
  cpp <- file.path(src_dir, paste0(short, ".cpp"))
  TMB::compile(cpp, ...)
  invisible(cpp)
}
make_tmb <- function(short, d = get_pdb_model(short)$data, compile = FALSE) {
  if (compile || !file.exists(file.path(src_dir, TMB::dynlib(short)))) compile_tmb(short)
  dyn.load(file.path(src_dir, TMB::dynlib(short)))
  tmb_data <- lapply(d, function(x) if (is.integer(x)) as.numeric(x) else x)
  TMB::MakeADFun(tmb_data, model_specs[[short]]$par(d), DLL = short, silent = TRUE)
}

## glmmTMB: returns the fitted model unless doFit = FALSE. With priors = TRUE,
## adds the model's glmmTMB_priors (ML becomes MAP).
make_glmmTMB <- function(short, d = get_pdb_model(short)$data, priors = FALSE, ...) {
  spec <- model_specs[[short]]
  glmmTMB::glmmTMB(spec$glmmTMB$formula, family = spec$glmmTMB$family,
                   data = spec$df(d),
                   priors = if (priors) spec$glmmTMB_priors, ...)
}

## brms: generated Stan code and data (compiled and run via cmdstanr/stanli)
make_brms_code <- function(short, d = get_pdb_model(short)$data) {
  spec <- model_specs[[short]]$brms
  df <- model_specs[[short]]$df(d)
  list(code = brms::stancode(spec$formula, data = df, family = spec$family,
                             prior = spec$prior),
       data = unclass(brms::standata(spec$formula, data = df, family = spec$family,
                                     prior = spec$prior)))
}
