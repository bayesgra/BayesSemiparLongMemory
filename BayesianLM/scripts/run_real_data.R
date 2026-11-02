file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript scripts/run_real_data.R")
repo_root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
source(file.path(repo_root, "scripts", "_common.R"))
opt <- cli_options()
check_options(opt, c("dataset", "runs", "burn", "thin", "H", "seed", "out", "taper"))
dataset <- text_option(opt, "dataset", "sp500")
if (!dataset %in% c("sp500", "sp500_covid", "temperature", "nile"))
  stop("dataset must be sp500, sp500_covid, temperature or nile")
dates <- NULL
if (dataset %in% c("sp500", "sp500_covid")) {
  raw <- read.csv(file.path(repo_root, "data", "SP500.csv"), stringsAsFactors = FALSE)
  raw$SP500 <- suppressWarnings(as.numeric(raw$SP500))
  raw$DATE <- as.Date(raw$DATE)
  raw <- raw[is.finite(raw$SP500) & !is.na(raw$DATE), ]
  raw <- raw[order(raw$DATE), ]
  if (anyDuplicated(raw$DATE) || any(raw$SP500 <= 0)) stop("Invalid price observations")
  series <- diff(log(raw$SP500))
  dates <- raw$DATE[-1L]
  if (dataset == "sp500_covid") {
    selected <- dates >= as.Date("2019-07-30") & dates <= as.Date("2020-10-31")
    series <- series[selected]
    dates <- dates[selected]
  }
  transformation <- "daily log-price differences; returns indexed by ending date"
} else if (dataset == "temperature") {
  raw <- read.csv(file.path(repo_root, "data", "HadCRUT.5.0.2.0.analysis.summary_series.global.monthly.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
  series <- as.numeric(raw[["Anomaly (deg C)"]])
  if (any(!is.finite(series))) stop("Missing temperature anomalies: no automatic gap deletion")
  dates <- as.Date(paste0(raw$Time, "-01"))
  if (anyNA(dates) || anyDuplicated(dates) || is.unsorted(dates)) stop("Invalid temperature dates")
  expected <- seq(dates[1L], by = "month", length.out = length(dates))
  if (!identical(dates, expected)) stop("Temperature record has gaps")
  series <- diff(series)
  dates <- dates[-1L]
  transformation <- "first differences of monthly HadCRUT anomalies, following supplied legacy code"
} else {
  if (!requireNamespace("longmemo", quietly = TRUE))
    stop("The original Nile data are longmemo::NileMin. Install longmemo separately; datasets::Nile is a different series.")
  env <- new.env()
  data("NileMin", package = "longmemo", envir = env)
  if (!exists("NileMin", env, inherits = FALSE)) stop("longmemo did not provide NileMin")
  series <- as.numeric(env$NileMin)
  transformation <- "untransformed longmemo::NileMin, following supplied legacy code"
}
runs <- number_option(opt, "runs", 30000L)
burn <- number_option(opt, "burn", min(6000L, floor(runs / 2)))
taper <- text_option(opt, "taper", "order1")
dat <- periodogram_data(series, taper = taper)
fit <- longmemory_ksb_sampler(dat$Y, dat$omega, H = number_option(opt, "H", 30L),
                             runs = runs, burn = burn, thin = number_option(opt, "thin", 5L),
                             seed = number_option(opt, "seed", 2026L), display = 0)
# Paper Section 6 reports an untrimmed GPH application; Section 5 uses first=5.
freq <- gph_estimate(series, first = 1L)
out <- text_option(opt, "out", file.path(repo_root, "results", dataset))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(fit = fit, periodogram = dat, gph = freq, series = series, dates = dates,
             transformation = transformation, options = opt, session = sessionInfo()),
        file.path(out, "analysis.rds"))
summary <- fit_summary(fit)
summary$method <- "BayesSemi"
summary <- rbind(summary, data.frame(estimate = freq$estimate, sd = freq$se,
                                     lower = freq$lower, upper = freq$upper, truth = NA_real_,
                                     covered = NA, method = "GPH_baseR"))
write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(data.frame(date = if (is.null(dates)) seq_along(series) else dates, value = series),
          file.path(out, "analyzed_series.csv"), row.names = FALSE)
save_diagnostic_plot(series, dat, fit, file.path(out, "diagnostics.pdf"), dataset)
cat("Saved", dataset, "analysis to", normalizePath(out), "\n")
cat("Fourier cutoff:", dat$cutoff, "; ordinates retained:", length(dat$Y), "\n")
cat("Check convergence and truncation sensitivity before interpreting intervals.\n")
