## Worker for R/bench_eval.R: setup and function/gradient timings for one
## engine on one model, run in a fresh R process. Returns a data frame with
## one row per metric (seconds; per call for fn/grad metrics).
##
## Engines: cmdstan, stanli, stanli_fast (posteriordb code); brms_cmdstan,
## brms_stanli (brms-generated code); rtmb, tmb; glmmTMB.
## The cmdstan and brms_cmdstan workers also make the evaluation points
## (results/points/), so they must run before the others.

source(here::here("R", "bench_utils.R"))

eval_worker <- function(engine, short, n_compile = 3, min_time = 0.1) {
  n_threads_1()
  set.seed(20261115)
  program <- if (grepl("^brms_", engine)) "brms" else "pdb"
  out <- function(metric, seconds, n_par = NA_integer_) {
    data.frame(model = short, program = program, engine = engine,
               metric = metric, seconds = seconds, n_par = n_par)
  }
  get_points <- function() readRDS(points_file(short, program))

  switch(engine,

    cmdstan = , brms_cmdstan = {
      p <- get_program(short, program)
      stanc <- file.path(cmdstanr::cmdstan_path(), "bin", "stanc")
      hpp <- tempfile(fileext = ".hpp")
      t_stanc <- time_reps(system2(stanc, c(p$stan_file, paste0("--o=", hpp)),
                                   stdout = FALSE), n_compile)
      t_compile <- time_reps(cmdstanr::cmdstan_model(p$stan_file, force_recompile = TRUE,
                                                     quiet = TRUE), n_compile)
      ## model methods (log_prob() etc.) are only built by a fresh compile
      t_methods <- system.time(
        mod <- cmdstanr::cmdstan_model(p$stan_file, compile_model_methods = TRUE,
                                       force_recompile = TRUE, quiet = TRUE)
      )[["elapsed"]]
      ## reference run -> evaluation points and gradient count
      fit <- mod$sample(data = p$data, chains = 1, iter_warmup = 1000,
                        iter_sampling = 1000, save_warmup = TRUE, seed = 1,
                        refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
      nlf <- as.numeric(posterior::extract_variable(
        fit$sampler_diagnostics(inc_warmup = TRUE), "n_leapfrog__"))
      idx <- select_points(nlf)
      pars <- names(mod$variables()$parameters)
      ## (unconstrain_draws() can't take a subset of iterations)
      da <- fit$draws(variables = pars, inc_warmup = TRUE)
      U <- t(sapply(idx, function(i) {
        r <- posterior::as_draws_rvars(posterior::subset_draws(da, iteration = i))
        fit$unconstrain_variables(lapply(r, function(x) drop(posterior::draws_of(x))))
      }))
      if (nrow(U) != n_points) U <- t(U)  ## sapply() simplifies 1-parameter models differently
      dir.create(points_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(list(U = U, iterations = idx, n_grad = sum(nlf)), points_file(short, program))
      rbind(
        if (program == "brms") out("codegen", p$gen_time),
        out("stanc", t_stanc), out("compile", t_compile),
        out("compile_methods", t_methods),
        out("n_grad_ref", sum(nlf), ncol(U)),
        out("fn", time_per_call(function(u) fit$log_prob(u, jacobian = TRUE), U, min_time), ncol(U)),
        out("grad", time_per_call(function(u) fit$grad_log_prob(u, jacobian = TRUE), U, min_time), ncol(U))
      )
    },

    stanli = , stanli_fast = , brms_stanli = {
      library(stanli)
      p <- get_program(short, program)
      fast <- engine == "stanli_fast"
      t_setup <- time_reps(stanli_model(file = p$stan_file, data = p$data, fast_math = fast),
                           max(6, n_compile))
      sl <- stanli_model(file = p$stan_file, data = p$data, fast_math = fast)
      U <- get_points()$U
      ## stanli's R interface has no log-density-only call
      rbind(out("setup", t_setup),
            out("fngrad", time_per_call(function(u) log_prob_grad(sl, u), U, min_time), ncol(U)))
    },

    rtmb = , tmb = {
      d <- get_pdb_model(short)$data
      if (engine == "tmb") {
        ## remove the compiled objects, so every compile starts from scratch
        rm_obj <- function() unlink(file.path(src_dir, paste0(short, c(".o", ".so"))))
        t_compile <- median(replicate(n_compile, {
          rm_obj()
          system.time(compile_tmb(short))[["elapsed"]]
        }))
        t_setup <- time_reps(make_tmb(short, d), max(6, n_compile))
        obj <- make_tmb(short, d)
      } else {
        t_compile <- NA_real_
        t_setup <- time_reps(make_rtmb(short, d), max(6, n_compile))
        obj <- make_rtmb(short, d)
      }
      U <- get_points()$U
      np <- ncol(U)
      rbind(
        if (engine == "tmb") out("compile", t_compile),
        out("setup", t_setup),
        out("fn", time_per_call(function(u) obj$fn(u), U, min_time), np),
        out("grad", time_per_call(function(u) obj$gr(u), U, min_time), np),
        ## obj$env$f is what tmbstan and StanEstimators call; it skips the R
        ## bookkeeping in obj$fn/obj$gr
        out("fn_raw", time_per_call(function(u) obj$env$f(u), U, min_time), np),
        out("grad_raw", time_per_call(function(u) obj$env$f(u, order = 1), U, min_time), np)
      )
    },

    glmmTMB = {
      d <- get_pdb_model(short)$data
      t_setup <- time_reps(make_glmmTMB(short, d, doFit = FALSE), max(6, n_compile))
      t_fit <- time_reps(make_glmmTMB(short, d), max(6, n_compile))
      fit <- make_glmmTMB(short, d)
      obj <- fit$obj
      ## Laplace objective (outer parameters, each call runs the inner
      ## optimization) at the MLE, and the joint density over all parameters
      ## (what tmbstan samples) at the mode
      par_outer <- matrix(fit$fit$par, nrow = 1)
      par_joint <- matrix(obj$env$last.par.best, nrow = 1)
      rbind(
        out("setup", t_setup), out("ml_fit", t_fit),
        out("laplace_fn", time_per_call(function(u) obj$fn(u), par_outer, min_time), ncol(par_outer)),
        out("laplace_grad", time_per_call(function(u) obj$gr(u), par_outer, min_time), ncol(par_outer)),
        out("fn_raw", time_per_call(function(u) obj$env$f(u), par_joint, min_time), ncol(par_joint)),
        out("grad_raw", time_per_call(function(u) obj$env$f(u, order = 1), par_joint, min_time),
            ncol(par_joint))
      )
    },
    stop("unknown engine: ", engine)
  )
}
