## Packages for the cmdstanr / stanli / brms / TMB / RTMB / glmmTMB benchmarks.
## Run once (needs internet access) from the project root:
##   Rscript get_pkgs.R
## Installs only what is missing, then installs CmdStan and the stanli runtime
## if they aren't already present.
cran <- "https://cloud.r-project.org"
options(repos = c(CRAN = cran))

## CRAN packages
pkgs <- c(
  "bench", "brms", "glmmTMB", "here", "posterior", "RhpcBLASctl",
  "RTMB", "TMB", "tibble", "tidyverse", "tmbstan"
)
pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs) > 0) install.packages(pkgs)

## packages from r-universe repositories (not on CRAN)
universe_pkgs <- c(
  cmdstanr       = "https://stan-dev.r-universe.dev",
  posteriordb    = "https://stan-dev.r-universe.dev",
  StanEstimators = "https://andrjohns.r-universe.dev",
  stanli         = "https://seantalts.r-universe.dev"
)
for (p in names(universe_pkgs)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = c(universe_pkgs[[p]], CRAN = cran))
  }
}

## CmdStan itself (used by cmdstanr and brms). stanli's published benchmarks
## compare against CmdStan 2.40.0, so require at least that version.
cmdstan_min <- "2.40.0"
cs_ver <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
if (is.null(cs_ver) || package_version(cs_ver) < cmdstan_min) {
  cmdstanr::check_cmdstan_toolchain(fix = TRUE)
  cmdstanr::install_cmdstan(version = cmdstan_min,
                            cores = max(1, parallel::detectCores() %/% 2, na.rm = TRUE))
}

## stanli downloads its prebuilt runtime on first use; fetch it now
stanli::stanli_install()

## check
for (p in c("glmmTMB", "TMB", "RTMB", "tmbstan", "StanEstimators",
            "cmdstanr", "brms", "stanli", "posteriordb", "posterior")) {
  cat(p, as.character(packageVersion(p)), "\n")
}
cat("CmdStan", cmdstanr::cmdstan_version(), "\n")
