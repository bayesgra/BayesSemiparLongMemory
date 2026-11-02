# Finite zero-start convolution approximation to stationary ARFIMA(1,d,1).
# (1-phi*B)(1-B)^d X_t = (1+theta*B) epsilon_t.
# A long warmup is discarded. This is NOT an exact stationary simulator;
# assess burn sensitivity, especially d close to 1/2.
simulate_arfima <- function(n, d, phi = 0, theta = 0,
                            innovations = c("gaussian", "student"), df = 5,
                            burn = max(5000L, 5L * n), seed = NULL) {
  innovations <- match.arg(innovations)
  scalar <- function(z) is.numeric(z) && length(z) == 1L && is.finite(z)
  if (!scalar(n) || n != floor(n) || n < 8 || !scalar(burn) ||
      burn < 0 || burn != floor(burn)) stop("n and burn must be valid integer counts")
  if (!scalar(d) || d <= -1 || d >= .5 || !scalar(phi) || abs(phi) >= 1 ||
      !scalar(theta) || abs(theta) >= 1) stop("invalid stationary/invertible parameters")
  if (!scalar(df) || df <= 2) stop("df must exceed two")
  if (!is.null(seed)) set.seed(seed)
  total <- as.integer(n + burn)
  eps <- if (innovations == "gaussian") rnorm(total) else rt(total, df) * sqrt((df - 2) / df)
  ma <- eps + theta * c(0, head(eps, -1L))
  short <- as.numeric(filter(ma, phi, method = "recursive", init = 0))
  coeff <- numeric(total)
  coeff[1L] <- 1
  for (k in seq_len(total - 1L)) coeff[k + 1L] <- coeff[k] * (k - 1 + d) / k
  fft_n <- 2^ceiling(log2(2 * total - 1))
  padded <- function(z) c(z, numeric(fft_n - length(z)))
  full <- Re(fft(fft(padded(short)) * fft(padded(coeff)), inverse = TRUE)) / fft_n
  ans <- full[burn + seq_len(n)]
  attr(ans, "simulation") <- list(n = n, d = d, phi = phi, theta = theta,
                                  innovations = innovations, df = df, burn = burn,
                                  method = "finite zero-start convolution with discarded warmup")
  ans
}
