file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript scripts/run_simulations.R")
repo_root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
source(file.path(repo_root, "scripts", "_common.R"))
opt <- cli_options()
check_options(opt, c("mode", "runs", "burn", "thin", "H", "seed", "reps", "out", "sim_burn"))
mode <- text_option(opt, "mode", "smoke")
if (!mode %in% c("smoke", "paper")) stop("mode must be smoke or paper")
paper <- mode == "paper"
runs <- number_option(opt, "runs", if (paper) 30000L else 200L)
burn <- number_option(opt, "burn", if (paper) 6000L else floor(runs / 2))
thin <- number_option(opt, "thin", if (paper) 5L else 1L)
H <- number_option(opt, "H", if (paper) 30L else 5L)
reps <- number_option(opt, "reps", if (paper) 200L else 2L)
seed <- number_option(opt, "seed", 6533L)
if (reps < 1L || reps != floor(reps)) stop("reps must be a positive integer")
out <- text_option(opt, "out", file.path(repo_root, "results", paste0("simulation_", mode)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "chains"), showWarnings = FALSE)
short <- data.frame(phi = c(0, .8, .5, .1, .1), theta = c(0, .1, .1, .5, .8))
grid <- if (paper) expand.grid(n = c(5000L, 10000L, 20000L),
                               d = c(.05, .1, .2, .3, .4, .49), scenario = seq_len(5L)) else
  data.frame(n = 512L, d = .2, scenario = 1L)
grid$phi <- short$phi[grid$scenario]
grid$theta <- short$theta[grid$scenario]
grid$innovations <- "gaussian"
heavy <- grid[grid$scenario == 1L, , drop = FALSE]
heavy$innovations <- "student"
grid <- rbind(grid, heavy)
write.csv(grid, file.path(out, "settings.csv"), row.names = FALSE)
saveRDS(list(options = opt, mode = mode, H = H, runs = runs, burn = burn, thin = thin,
             reps = reps, seed = seed, session = sessionInfo(),
             caveat = "Finite zero-start ARFIMA simulation; inspect warmup sensitivity and MCMC convergence"),
        file.path(out, "metadata.rds"))
results <- list()
idx <- 0L
for (g in seq_len(nrow(grid))) for (r in seq_len(reps)) {
  setting <- grid[g, ]
  simulation_seed <- seed + 2L * ((g - 1L) * reps + r)
  series <- simulate_arfima(setting$n, setting$d, setting$phi, setting$theta,
                            setting$innovations, seed = simulation_seed,
                            burn = number_option(opt, "sim_burn", max(5000L, 5L * setting$n)))
  dat <- periodogram_data(series)
  fit <- longmemory_ksb_sampler(dat$Y, dat$omega, H = H, runs = runs, burn = burn,
                               thin = thin, seed = simulation_seed + 1L, display = 0)
  semi <- fit_summary(fit, setting$d)
  semi$method <- "BayesSemi"
  freq <- gph_estimate(series, first = 5L)
  frequentist <- data.frame(estimate = freq$estimate, sd = freq$se, lower = freq$lower,
                            upper = freq$upper, truth = setting$d,
                            covered = freq$lower <= setting$d && setting$d <= freq$upper,
                            method = "GPH_baseR")
  idx <- idx + 1L
  results[[idx]] <- cbind(setting[rep(1L, 2L), , drop = FALSE], replication = r, seed = simulation_seed,
                          rbind(semi, frequentist))
  # Checkpoint completed summaries after every fit; rerun to a fresh directory.
  write.csv(do.call(rbind, results), file.path(out, "estimates.csv"), row.names = FALSE)
  saveRDS(list(fit = fit, setting = setting, replication = r, periodogram = dat),
          file.path(out, "latest_fit.rds"))
  saveRDS(list(d = fit$d, c = fit$c, diagnostics = fit$diagnostics,
               acceptance = fit$diagnostics$acceptance, setting = setting, replication = r,
               simulation_seed = simulation_seed, sampler_seed = simulation_seed + 1L),
          file.path(out, "chains", sprintf("setting_%03d_rep_%03d.rds", g, r)))
  cat("Completed setting", g, "of", nrow(grid), "replication", r, "of", reps, "\n")
}
cat("No BayesParam results are generated: its legacy implementation requires separate validation.\n")
cat("Paper mode schedules new runs; it does not reproduce or certify the manuscript's numerical tables.\n")
