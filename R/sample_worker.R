## Worker for R/bench_sample.R: one sampling run (4 chains run sequentially,
## 1000 warmup + 1000 draws) for one engine on one model, in a fresh R
## process. Returns a one-row data frame of times, ESS and diagnostics.
##
## Engines: cmdstan, stanli, stanli_fast (posteriordb code); brms_cmdstan,
## brms_stanli (brms-generated code); brms_brm (brm() itself, cmdstanr
## backend); tmbstan_rtmb, tmbstan_tmb, tmbstan_glmmTMB,
## tmbstan_glmmTMB_prior (glmmTMB with priors); stanest_rtmb (StanEstimators
## with the RTMB objective).
##
## Chains run one after another everywhere, so the summed chain time is the
## CPU cost and no engine benefits from parallel chains. ESS is computed for
## the sampled parameters (Stan's parameters block, or all TMB parameters);
## bulk ESS is rank-based, so it doesn't depend on the constrained vs
## unconstrained scale.

source(here::here("R", "bench_utils.R"))

sample_worker <- function(engine, short, seed, chains = 4, warmup = 1000, samples = 1000) {
  n_threads_1()
  set.seed(seed)
  program <- if (grepl("^brms_", engine)) "brms" else "pdb"

  ## each branch sets: wall, t_warmup, t_sampling (summed over chains),
  ## d (draws of sampled parameters, no warmup), n_leapfrog (warmup +
  ## sampling, summed), n_div (sampling)
  wall <- NA_real_
  t_compile <- NA_real_  ## only brms_brm compiles inside the sampling run
  switch(engine,

    cmdstan = , brms_cmdstan = {
      p <- get_program(short, program)
      mod <- cmdstanr::cmdstan_model(p$stan_file, quiet = TRUE)
      wall <- system.time(
        fit <- mod$sample(data = p$data, chains = chains, parallel_chains = 1,
                          iter_warmup = warmup, iter_sampling = samples,
                          save_warmup = TRUE, seed = seed, refresh = 0,
                          show_messages = FALSE, show_exceptions = FALSE)
      )[["elapsed"]]
      tm <- fit$time()$chains
      t_warmup <- sum(tm$warmup); t_sampling <- sum(tm$sampling)
      d <- fit$draws(variables = names(mod$variables()$parameters))
      sd <- fit$sampler_diagnostics(inc_warmup = TRUE)
      n_leapfrog <- sum(posterior::extract_variable(sd, "n_leapfrog__"))
      n_div <- sum(posterior::extract_variable(fit$sampler_diagnostics(), "divergent__"))
    },

    stanli = , stanli_fast = , brms_stanli = {
      library(stanli)
      p <- get_program(short, program)
      sl <- stanli_model(file = p$stan_file, data = p$data,
                         fast_math = engine == "stanli_fast")
      wall <- system.time(
        fit <- sample_model(sl, chains = chains, parallel_chains = 1, seed = seed,
                            warmup = warmup, samples = samples, save_warmup = TRUE,
                            refresh = 0)
      )[["elapsed"]]
      tm <- get_elapsed_time(fit)
      t_warmup <- sum(tm[, "warmup"]); t_sampling <- sum(tm[, "sample"])
      d <- select_pars(as_draws_array(fit), stan_param_names(p$stan_file))
      sp <- get_sampler_params(fit, inc_warmup = TRUE)
      n_leapfrog <- sum(sapply(sp, function(x) sum(x[, "n_leapfrog__"])))
      n_div <- sum(sapply(sp, function(x) sum(x[-seq_len(warmup), "divergent__"])))
    },

    tmbstan_rtmb = , tmbstan_tmb = , tmbstan_glmmTMB = , tmbstan_glmmTMB_prior = {
      d0 <- get_pdb_model(short)$data
      obj <- switch(engine,
                    tmbstan_rtmb = make_rtmb(short, d0),
                    tmbstan_tmb = make_tmb(short, d0),
                    tmbstan_glmmTMB = make_glmmTMB(short, d0)$obj,
                    tmbstan_glmmTMB_prior = make_glmmTMB(short, d0, priors = TRUE)$obj)
      wall <- system.time(
        fit <- tmbstan::tmbstan(obj, chains = chains, cores = 1, iter = warmup + samples,
                                warmup = warmup, seed = seed, refresh = 0,
                                init = "random")
      )[["elapsed"]]
      tm <- rstan::get_elapsed_time(fit)
      t_warmup <- sum(tm[, "warmup"]); t_sampling <- sum(tm[, "sample"])
      d <- posterior::as_draws_array(fit)
      d <- posterior::subset_draws(d, variable = setdiff(posterior::variables(d), "lp__"))
      sp <- rstan::get_sampler_params(fit, inc_warmup = TRUE)
      n_leapfrog <- sum(sapply(sp, function(x) sum(x[, "n_leapfrog__"])))
      n_div <- sum(sapply(sp, function(x) sum(x[-seq_len(warmup), "divergent__"])))
    },

    brms_brm = {
      ## the brms user's path: brm() does its own data processing and
      ## compilation (chains = 0 compiles without sampling), then samples from
      ## the compiled model and builds the brmsfit
      spec <- model_specs[[short]]$brms
      df <- model_specs[[short]]$df(get_pdb_model(short)$data)
      t_compile <- system.time(
        fit0 <- brms::brm(spec$formula, data = df, family = spec$family, prior = spec$prior,
                          backend = "cmdstanr", chains = 0, silent = 2, refresh = 0)
      )[["elapsed"]]
      wall <- system.time(
        fit <- stats::update(fit0, recompile = FALSE, chains = chains, cores = 1,
                             iter = warmup + samples, warmup = warmup, seed = seed,
                             silent = 2, refresh = 0)
      )[["elapsed"]]
      tm <- rstan::get_elapsed_time(fit$fit)
      t_warmup <- sum(tm[, "warmup"]); t_sampling <- sum(tm[, "sample"])
      d <- posterior::as_draws_array(fit)
      d <- posterior::subset_draws(d, variable = setdiff(posterior::variables(d),
                                                         c("lp__", "lprior")))
      np <- brms::nuts_params(fit)
      ## brm() doesn't keep warmup diagnostics, so the warmup gradient count
      ## (and hence time_per_grad) is unavailable
      n_leapfrog <- NA_real_
      n_div <- sum(np$Value[np$Parameter == "divergent__"])
    },

    stanest_rtmb = {
      obj <- make_rtmb(short, get_pdb_model(short)$data)
      ## obj$env$f skips the R bookkeeping of obj$fn/obj$gr (as tmbstan does)
      wall <- system.time(
        fit <- StanEstimators::stan_sample(
          fn = function(x) -obj$env$f(x),
          grad_fun = function(x) -as.numeric(obj$env$f(x, order = 1)),
          par_inits = obj$par, num_chains = chains, parallel_chains = 1,
          num_warmup = warmup, num_samples = samples, save_warmup = TRUE,
          seed = seed, refresh = 0, quiet = TRUE)
      )[["elapsed"]]
      t_warmup <- sum(sapply(fit@timing, `[[`, "warmup"))
      t_sampling <- sum(sapply(fit@timing, `[[`, "sampling"))
      keep <- function(x) posterior::subset_draws(x, iteration = warmup + seq_len(samples))
      d <- posterior::as_draws_array(fit@draws)
      d <- keep(posterior::subset_draws(d, variable = setdiff(posterior::variables(d), "lp__")))
      dg <- posterior::as_draws_array(fit@diagnostics)
      n_leapfrog <- sum(posterior::extract_variable(dg, "n_leapfrog__"))
      n_div <- sum(posterior::extract_variable(keep(dg), "divergent__"))
    },
    stop("unknown engine: ", engine)
  )

  stopifnot(posterior::niterations(d) == samples, posterior::nchains(d) == chains)
  s <- summarize_draws_ess(d)
  t_total <- t_warmup + t_sampling
  data.frame(model = short, program = program, engine = engine, seed = seed,
             t_compile = t_compile, wall = wall, t_warmup = t_warmup, t_sampling = t_sampling, t_total = t_total,
             s, n_leapfrog = n_leapfrog, n_div = n_div,
             time_per_grad = t_total / n_leapfrog,
             time_per_1000_ess = 1000 * t_total / s$ess_bulk_min,
             sampling_time_per_1000_ess = 1000 * t_sampling / s$ess_bulk_min)
}
