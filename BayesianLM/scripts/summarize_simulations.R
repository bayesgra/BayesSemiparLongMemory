file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run with Rscript scripts/summarize_simulations.R")
repo_root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg))))
source(file.path(repo_root, "scripts", "_common.R"))
opt <- cli_options()
check_options(opt, c("input", "out"))
input <- text_option(opt, "input", file.path(repo_root, "results", "simulation_smoke", "estimates.csv"))
out <- text_option(opt, "out", dirname(input))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
est <- read.csv(input, stringsAsFactors = FALSE)
required <- c("n", "d", "phi", "theta", "innovations", "method", "estimate", "sd", "covered")
if (!all(required %in% names(est))) stop("Input is not an estimates.csv produced by this project")
keys <- c("n", "d", "phi", "theta", "innovations", "method")
groups <- split(est, interaction(est[keys], drop = TRUE, lex.order = TRUE))
summary <- do.call(rbind, lapply(groups, function(z) {
  cbind(z[1L, keys, drop = FALSE], replications = nrow(z),
        bias = mean(z$estimate - z$d), rmse = sqrt(mean((z$estimate - z$d)^2)),
        coverage = mean(z$covered), mean_sd = mean(z$sd))
}))
rownames(summary) <- NULL
write.csv(summary, file.path(out, "performance.csv"), row.names = FALSE)
pdf(file.path(out, "boxplots.pdf"), width = 10, height = 7)
for (n in unique(est$n)) for (innovation in unique(est$innovations)) {
  z <- est[est$n == n & est$innovations == innovation, ]
  if (!nrow(z)) next
  scenarios <- unique(z[c("phi", "theta")])
  par(mfrow = c(ceiling(nrow(scenarios) / 2), 2), mar = c(5, 4, 3, 1))
  for (s in seq_len(nrow(scenarios))) {
    zz <- z[z$phi == scenarios$phi[s] & z$theta == scenarios$theta[s], ]
    boxplot(estimate ~ interaction(d, method, sep = ":", lex.order = TRUE), data = zz,
            las = 2, cex.axis = .65, xlab = "d : method", ylab = "Estimate of d",
            main = paste("n =", n, innovation, "phi =", scenarios$phi[s], "theta =", scenarios$theta[s]))
    abline(h = sort(unique(zz$d)), col = "grey60", lty = 2)
  }
}
dev.off()
cat("Saved performance.csv and boxplots.pdf to", normalizePath(out), "\n")
