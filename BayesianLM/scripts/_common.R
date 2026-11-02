# Sourced only by entry-point scripts after they resolve repository root.
source(file.path(repo_root, "R", "periodogram.R"))
source(file.path(repo_root, "R", "simulate_arfima.R"))
source(file.path(repo_root, "R", "longmemory_ksb.R"))
cli_options <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  if (any(!grepl("^--[A-Za-z_]+=[^=]+$", a)))
    stop("Use --name=value arguments (see README)")
  setNames(as.list(sub("^[^=]+=", "", a)), sub("=.+$", "", sub("^--", "", a)))
}
number_option <- function(options, name, default) {
  if (is.null(options[[name]])) return(default)
  z <- suppressWarnings(as.numeric(options[[name]]))
  if (length(z) != 1L || !is.finite(z)) stop(paste("Invalid option", name))
  z
}
text_option <- function(options, name, default) {
  if (is.null(options[[name]])) default else options[[name]]
}
check_options <- function(options, allowed) {
  bad <- setdiff(names(options), allowed)
  if (length(bad)) stop(paste("Unknown option(s):", paste(bad, collapse = ", ")))
}
fit_summary <- function(fit, truth = NA_real_) {
  interval <- quantile(fit$d, c(.025, .975), names = FALSE)
  data.frame(estimate = mean(fit$d), sd = sd(fit$d), lower = interval[1],
             upper = interval[2], truth = truth,
             covered = if (is.na(truth)) NA else interval[1] <= truth && truth <= interval[2])
}
save_diagnostic_plot <- function(series, dat, fit, filename, title = "Analysis") {
  pdf(filename, width = 10, height = 8)
  on.exit(dev.off())
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  plot(series, type = "l", xlab = "Observation", ylab = "Value", main = title)
  acf(series, lag.max = min(40L, length(series) - 1L), main = "Sample autocorrelation")
  plot(dat$x, dat$Y, pch = 16, cex = .5, xlab = "log(4 sin²(omega/2))",
       ylab = "Log periodogram", main = "Retained Fourier ordinates")
  abline(a = mean(fit$c), b = -mean(fit$d), col = "steelblue", lwd = 2)
  hist(fit$d, probability = TRUE, breaks = "FD", col = "grey85", border = "white",
       main = "Marginal pseudo-posterior draws", xlab = "d")
  abline(v = mean(fit$d), col = "steelblue", lwd = 2)
  abline(v = quantile(fit$d, c(.025, .975)), lty = 2)
}
