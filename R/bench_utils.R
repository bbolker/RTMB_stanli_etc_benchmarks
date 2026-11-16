## Helpers and per-process workers for the timing scripts.
##
## Each (engine, model) measurement runs in a fresh R process (callr::r): the
## compiled model methods of cmdstanr segfault when loaded into a session that
## already has rstan/tmbstan or StanEstimators (each brings its own Stan Math
## and TBB), and fresh processes also keep one engine's state from affecting
## another's timings.

source(here::here("R", "models.R"))

points_dir <- here::here("results", "points")
n_points <- 10

## ---- programs: posteriordb Stan code, or brms-generated code ----

## Stan file and data for "pdb" (posteriordb) or "brms" (brms-generated) code
get_program <- function(short, program = c("pdb", "brms")) {
  program <- match.arg(program)
  if (program == "pdb") {
    m <- get_pdb_model(short)
    return(list(stan_file = m$stan_file, data = m$data, gen_time = NA_real_))
  }
  gen_time <- system.time(bc <- make_brms_code(short))[["elapsed"]]
  stan_file <- file.path(model_dir, paste0(short, "_brms.stan"))
  ## rewrite only on change, so cmdstanr can reuse an up-to-date executable
  if (!file.exists(stan_file) || !identical(readLines(stan_file), strsplit(bc$code, "\n")[[1]])) {
    writeLines(bc$code, stan_file)
  }
  list(stan_file = stan_file, data = bc$data, gen_time = gen_time)
}

## names of the Stan parameters-block variables (without indices)
stan_param_names <- function(stan_file) {
  names(cmdstanr::cmdstan_model(stan_file, compile = FALSE)$variables()$parameters)
}

## ---- timing ----

## Mean over the rows of U (points) of the median time per call of f(u).
## As in stanli's protocol: evaluate once untimed, warm up for >= 20 ms, then
## time for >= min_time per point (bench::mark calibrates the batching).
time_per_call <- function(f, U, min_time = 0.1) {
  per_point <- apply(U, 1, function(u) {
    f(u)
    t0 <- proc.time()[["elapsed"]]
    while (proc.time()[["elapsed"]] - t0 < 0.02) f(u)
    b <- bench::mark(f(u), min_time = min_time, max_iterations = 1e7,
                     check = FALSE, filter_gc = FALSE, time_unit = "s")
    as.numeric(b$median)
  })
  mean(per_point)
}

## median elapsed time over `reps` evaluations of `expr`
time_reps <- function(expr, reps) {
  expr <- substitute(expr)
  env <- parent.frame()
  median(replicate(reps, system.time(eval(expr, env))[["elapsed"]]))
}

## ---- evaluation points (stanli's protocol) ----

## Pick n points from a 1-chain reference run with saved warmup: iteration i
## is point k (0..n-1) when it is the first whose cumulative gradient count
## reaches (k + 0.5)/n of the total. Returns the iteration indices.
select_points <- function(n_leapfrog, n = n_points) {
  cum <- cumsum(n_leapfrog)
  sapply((seq_len(n) - 0.5) / n, function(q) which(cum >= q * cum[length(cum)])[1])
}

points_file <- function(short, program) {
  file.path(points_dir, sprintf("%s_%s.rds", short, program))
}

## ---- sampling summaries ----

## draws_array/draws_df of the sampled parameters (no warmup) -> ESS etc.
summarize_draws_ess <- function(d) {
  s <- posterior::summarise_draws(d, "rhat", "ess_bulk", "ess_tail")
  data.frame(n_par = nrow(s),
             ess_bulk_min = min(s$ess_bulk), ess_bulk_med = median(s$ess_bulk),
             ess_tail_min = min(s$ess_tail), rhat_max = max(s$rhat))
}

## keep draws of variables whose base name is in `pars`
select_pars <- function(d, pars) {
  v <- posterior::variables(d)
  posterior::subset_draws(d, variable = v[sub("\\[.*", "", v) %in% pars])
}

n_threads_1 <- function() {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}
