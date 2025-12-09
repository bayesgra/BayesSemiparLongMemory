d0_vals <- c(0.05, 0.10, 0.15, 0.20,
             0.25, 0.30, 0.35, 0.40,
             0.45, 0.49)
K <- length(d0_vals)
R <- 50   # replications

d_gph_mat  <- matrix(NA, nrow = K, ncol = R)
d_semi_mat <- matrix(NA, nrow = K, ncol = R)
d_filt_mat <- matrix(NA, nrow = K, ncol = R)

for (k in 1:K) {
  file <- sprintf("BayesARFIMA-%03d.RData", round(100*d0_vals[k]))
  load(file)
  
  d_gph_mat[k, ]  <- d_gph_vec
  d_semi_mat[k, ] <- d_semi_vec
  d_filt_mat[k, ] <- d_filtered
}

plot_d_boxplots <- function(d0_vals, d_gph_mat, d_semi_mat, d_filt_mat) {
  
  K <- length(d0_vals)
  R <- ncol(d_gph_mat)
  
  cols <- c("blue", "red", "cyan")   # Freq, Semi, Filtered
  
  # X-axis positions
  x_base <- 1:K
  offset <- c(-0.25, 0, 0.25)
  
  # Empty plot
  plot(NA, xlim = c(0.5, K + 0.5),
       ylim = range(d_gph_mat, d_semi_mat, d_filt_mat),
       xaxt = "n", xlab = expression(d[0]), ylab = "Estimate")
  
  axis(1, at = x_base, labels = d0_vals)
  
  # Loop over each true d0
  for (k in 1:K) {
    # all estimates for the kth true d
    y1 <- d_semi_mat[k, ]
    y2 <- d_gph_mat[k, ]
    y3 <- d_filt_mat[k, ]
    
    # positions for the three boxplots
    x_pos <- x_base[k] + offset
    
    boxplot(y1, at = x_pos[1], col = "red",   add = TRUE, outline = TRUE)
    boxplot(y2, at = x_pos[2], col = "blue",  add = TRUE, outline = TRUE)
    boxplot(y3, at = x_pos[3], col = "cyan",  add = TRUE, outline = TRUE)
  }
  
#  legend("topleft",
#         legend = c("Frequentist", "Semiparametric", "Filtered Bayesian"),
#         fill   = cols, bty = "n")
}

plot_d_boxplots(d0_vals, d_gph_mat, d_semi_mat, d_filt_mat)
