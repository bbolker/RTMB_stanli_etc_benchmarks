## Sampling timings and ESS for every engine x model x replicate, each run in
## a fresh R process (see R/bench_utils.R). Writes results/sample_timings.rds.
## Compiled models are reused when up to date (run R/bench_eval.R first, or
## the first replicate includes compilation in its wall time).
##
## Environment variables: MODELS, ENGINES (comma-separated; default all),
## N_REP (replicates per engine x model, with different seeds; default 3).
source(here::here("R", "models.R"))

env_list <- function(var, default) {
  x <- Sys.getenv(var)
  if (nzchar(x)) strsplit(x, ",")[[1]] else default
}
models <- env_list("MODELS", model_names)
all_engines <- c("cmdstan", "stanli", "stanli_fast", "brms_cmdstan", "brms_stanli",
                 "brms_brm", "tmbstan_rtmb", "tmbstan_tmb", "tmbstan_glmmTMB",
                 "tmbstan_glmmTMB_prior", "stanest_rtmb")
engines <- intersect(all_engines, env_list("ENGINES", all_engines))
n_rep <- as.integer(Sys.getenv("N_REP", "3"))
seeds <- 20261115 + seq_len(n_rep)

run_one <- function(engine, short, seed) {
  cat(format(Sys.time(), "%H:%M:%S"), short, engine, seed, "\n")
  tryCatch(
    callr::r(function(engine, short, seed) {
      source(here::here("R", "sample_worker.R"))
      suppressWarnings(sample_worker(engine, short, seed))
    }, args = list(engine = engine, short = short, seed = seed),
    env = c(callr::rcmd_safe_env(), OPENBLAS_NUM_THREADS = "1",
            OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")),
    error = function(e) {
      message("FAILED ", short, " ", engine, " ", seed, ": ", conditionMessage(e))
      NULL
    })
}

## replicates outermost, so the engines alternate in time rather than each
## running all its replicates back to back
grid <- expand.grid(engine = engines, short = models, seed = seeds,
                    stringsAsFactors = FALSE)
res <- do.call(rbind, Map(run_one, grid$engine, grid$short, grid$seed))
rownames(res) <- NULL
res$machine <- Sys.info()[["nodename"]]

out_file <- here::here("results", "sample_timings.rds")
if (file.exists(out_file)) {
  old <- readRDS(out_file)
  key <- function(x) paste(x$model, x$engine, x$seed)
  res <- rbind(old[!key(old) %in% key(res), ], res)
}
saveRDS(res, out_file)
print(res[, c("model", "engine", "seed", "t_total", "ess_bulk_min", "n_div",
              "time_per_grad", "time_per_1000_ess")], digits = 3)
