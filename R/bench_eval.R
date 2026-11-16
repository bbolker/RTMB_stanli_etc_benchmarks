## Setup, function and gradient timings for every engine x model, each in a
## fresh R process (see R/bench_utils.R). Writes results/eval_timings.rds.
##
## Environment variables: MODELS (comma-separated short names; default all),
## ENGINES (default all), N_COMPILE (repeated compiles per program; default 3).
## With the defaults this takes roughly 30-60 min, mostly CmdStan compiles.
source(here::here("R", "models.R"))

env_list <- function(var, default) {
  x <- Sys.getenv(var)
  if (nzchar(x)) strsplit(x, ",")[[1]] else default
}
models <- env_list("MODELS", model_names)
## the cmdstan engines first: they write the evaluation points the others use
all_engines <- c("cmdstan", "brms_cmdstan", "stanli", "stanli_fast", "brms_stanli",
                 "rtmb", "tmb", "glmmTMB")
engines <- intersect(all_engines, env_list("ENGINES", all_engines))
n_compile <- as.integer(Sys.getenv("N_COMPILE", "3"))

run_one <- function(engine, short) {
  cat(format(Sys.time(), "%H:%M:%S"), short, engine, "\n")
  tryCatch(
    callr::r(function(engine, short, n_compile) {
      source(here::here("R", "eval_worker.R"))
      eval_worker(engine, short, n_compile)
    }, args = list(engine = engine, short = short, n_compile = n_compile),
    env = c(callr::rcmd_safe_env(), OPENBLAS_NUM_THREADS = "1",
            OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")),
    error = function(e) {
      message("FAILED ", short, " ", engine, ": ", conditionMessage(e))
      NULL
    })
}

res <- do.call(rbind, lapply(models, function(short) {
  do.call(rbind, lapply(engines, run_one, short = short))
}))

res$machine <- Sys.info()[["nodename"]]
out_file <- here::here("results", "eval_timings.rds")
## merge with earlier results, so partial runs (MODELS/ENGINES) update them
if (file.exists(out_file)) {
  old <- readRDS(out_file)
  old <- old[!(paste(old$model, old$engine) %in% paste(res$model, res$engine)), ]
  res <- rbind(old, res)
}
saveRDS(res, out_file)
print(res, digits = 3)
