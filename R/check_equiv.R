## Check that cmdstanr, stanli, RTMB and TMB compute the same log density and
## gradient (on Stan's unconstrained scale) at a set of points from a short
## cmdstanr run. Log densities may differ by a constant (Stan's `~` drops
## normalizing constants), so the lp check uses deviations from the mean
## difference across points. Gradients must match directly.
source(here::here("R", "models.R"))
library(cmdstanr)
library(stanli)
set.seed(20261115)

n_points <- 10

check_one <- function(short) {
  m <- get_pdb_model(short)
  ## cmdstanr only builds the model methods (log_prob() etc.) during a fresh
  ## compile; with an up-to-date executable they fail, hence force_recompile
  mod <- cmdstan_model(m$stan_file, compile_model_methods = TRUE,
                       force_recompile = TRUE, quiet = TRUE)
  fit <- mod$sample(data = m$data, chains = 1, iter_warmup = 500,
                    iter_sampling = 10 * n_points, thin = 10, seed = 20261115,
                    refresh = 0, show_messages = FALSE)
  ## one unconstrained point per row
  U <- unclass(fit$unconstrain_draws(format = "draws_matrix"))
  U <- matrix(U, nrow = nrow(U))
  sl <- stanli_model(file = m$stan_file, data = m$data)
  rt <- make_rtmb(short, m$data)
  tm <- make_tmb(short, m$data)
  stopifnot(ncol(U) == length(rt$par), ncol(U) == length(tm$par))

  res <- lapply(seq_len(nrow(U)), function(i) {
    u <- U[i, ]
    g_stan <- fit$grad_log_prob(u, jacobian = TRUE)
    s <- log_prob_grad(sl, u)
    list(lp = c(stan = fit$log_prob(u, jacobian = TRUE), stanli = s$lp,
                rtmb = -rt$fn(u), tmb = -tm$fn(u)),
         grad = cbind(stan = as.numeric(g_stan), stanli = s$grad,
                      rtmb = -as.numeric(rt$gr(u)), tmb = -as.numeric(tm$gr(u))))
  })
  lp <- do.call(rbind, lapply(res, `[[`, "lp"))
  engines <- c("stanli", "rtmb", "tmb")
  data.frame(
    model = short, n_par = ncol(U), engine = engines,
    lp_offset = sapply(engines, \(e) mean(lp[, e] - lp[, "stan"])),
    ## lp after removing the constant offset
    lp_maxdev = sapply(engines, \(e) {
      off <- mean(lp[, e] - lp[, "stan"])
      max(scaled_err(lp[, e] - off, lp[, "stan"]))
    }),
    grad_maxrel = sapply(engines, \(e) max(sapply(res, \(r) {
      max(scaled_err(r$grad[, e], r$grad[, "stan"]))
    }))),
    row.names = NULL
  )
}

## stanli's numerical gate: abs(a-b) / max(abs(a), abs(b), 1), per component
scaled_err <- function(a, b) abs(a - b) / pmax(abs(a), abs(b), 1)

equiv <- do.call(rbind, lapply(model_names, check_one))
print(equiv, digits = 3)
tol <- 1e-9
equiv$ok <- equiv$lp_maxdev < tol & equiv$grad_maxrel < tol
dir.create(here::here("results"), showWarnings = FALSE)
saveRDS(equiv, here::here("results", "equiv_check.rds"))
if (!all(equiv$ok)) stop("log density / gradient mismatch: see rows with ok = FALSE")
cat("all engines agree (tolerance", tol, ")\n")
