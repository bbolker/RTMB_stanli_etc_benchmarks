ALL: README.md

R = Rscript

## single-threaded linear algebra everywhere, so timings measure the engines
## rather than BLAS/OpenMP threading (exported to every Rscript call)
export OPENBLAS_NUM_THREADS = 1
export OMP_NUM_THREADS = 1
export MKL_NUM_THREADS = 1

MODEL_R = R/models.R $(wildcard src/*.cpp src/*.hpp)

.PHONY: get_pkgs
get_pkgs:
	$(R) get_pkgs.R

## log density and gradient agreement across cmdstanr, stanli, RTMB and TMB
## (compiles every Stan and TMB model; a few minutes)
results/equiv_check.rds: R/check_equiv.R $(MODEL_R)
	$(R) $<

BENCH_R = R/bench_utils.R $(MODEL_R)

## setup, function and gradient timings (about 30-60 min with N_COMPILE=3,
## mostly CmdStan and TMB compiles)
results/eval_timings.rds: R/bench_eval.R R/eval_worker.R $(BENCH_R)
	$(R) $<

## sampling timings and ESS (N_REP replicates, default 3); uses the compiled
## models from the evaluation timings
results/sample_timings.rds: R/bench_sample.R R/sample_worker.R $(BENCH_R) results/eval_timings.rds
	$(R) $<

.PHONY: bench
bench: results/eval_timings.rds results/sample_timings.rds

README.md: README.qmd results/equiv_check.rds
	quarto render $<
