# Optional conditional Gaussian ARFIMA benchmark, implemented using base R.
# This is NOT an implementation of the exact filtered method of Gabor & Grazian.
# Model: (1 - phi B)(1 - B)^d Y_t = (1 + theta B) epsilon_t.
# All presample values are set to zero. This finite-past conditional likelihood
# approximates the stationary likelihood and can be sensitive to initialization.

arfima_fractional_difference <- function(y, d) {
  n <- length(y)
  if (n == 1L) return(as.numeric(y))
  k <- seq_len(n - 1L)
  # Coefficients of (1-B)^d; cumprod correctly handles d=0.
  weights <- c(1, cumprod((k - 1 - d) / k))
  size <- 2^ceiling(log2(2 * n - 1))
  a <- c(as.numeric(y), rep(0, size - n))
  b <- c(weights, rep(0, size - n))
  Re(fft(fft(a) * fft(b), inverse = TRUE))[seq_len(n)] / size
}

arfima_conditional_innovations <- function(y, d, phi = 0, theta = 0) {
  u <- arfima_fractional_difference(y, d)
  z <- u - phi * c(0, head(u, -1L))
  # (1 + theta B) epsilon = z, hence epsilon_t=z_t-theta*epsilon_{t-1}.
  if (theta == 0) return(z)
  as.numeric(stats::filter(z, filter = -theta, method = "recursive", init = 0))
}

run_filtered_mcmc_arfima <- function(
    Y, p = 0, q = 0, n_iter = 30000, burn = 6000, thin = 5,
    d_sd = 0.025, arma_sd = 0.05, alpha = 0.1, beta = 0.1,
    seed = NULL, verbose = FALSE, init = c(d = 0, phi = 0, theta = 0)) {
  # Legacy function name retained only to simplify migration. Read method below.
  integer_scalar <- function(x) length(x) == 1L && is.numeric(x) &&
    is.finite(x) && x == floor(x)
  positive_scalar <- function(x) length(x) == 1L && is.numeric(x) &&
    is.finite(x) && x > 0
  if (!is.numeric(Y) || !is.null(dim(Y)) || length(Y) < 8L || any(!is.finite(Y)))
    stop("Y must be a finite numeric vector of at least 8 observations.")
  if (!integer_scalar(p) || !integer_scalar(q) || !p %in% 0:1 || !q %in% 0:1)
    stop("This benchmark supports known orders p,q in {0,1} only.")
  if (!integer_scalar(n_iter) || !integer_scalar(burn) || !integer_scalar(thin) ||
      n_iter < 1 || burn < 0 || burn >= n_iter || thin < 1 ||
      floor((n_iter - burn) / thin) < 1)
    stop("Require integer n_iter > burn >= 0, thin >= 1 and at least one saved draw.")
  if (!all(vapply(list(d_sd, arma_sd, alpha, beta), positive_scalar, logical(1))))
    stop("Proposal standard deviations and inverse-gamma parameters must be positive.")
  if (!is.numeric(init) || length(init) != 3L || any(!is.finite(init)))
    stop("init must contain finite d, phi, theta, in that order.")
  state <- setNames(as.numeric(init), c("d", "phi", "theta"))
  if (p == 0) state["phi"] <- 0
  if (q == 0) state["theta"] <- 0
  bounds <- c(d = 0.5, phi = 1, theta = 1)
  if (any(abs(state) >= bounds)) stop("Initial d, phi or theta is outside its support.")
  if (!is.null(seed)) {
    if (!integer_scalar(seed) || seed < 0 || seed > .Machine$integer.max)
      stop("seed must be NULL or a nonnegative integer.")
    set.seed(seed)
  }
  n <- length(Y)
  active <- c("d", if (p == 1) "phi", if (q == 1) "theta")
  proposal_sd <- c(d = d_sd, phi = arma_sd, theta = arma_sd)
  # The filtering map is lower triangular with unit diagonal, so its Jacobian
  # determinant is one. Integrating sigma^2 ~ IG(alpha,beta) gives this target.
  evaluate <- function(par) {
    residual <- arfima_conditional_innovations(Y, par["d"], par["phi"], par["theta"])
    sse <- sum(residual^2)
    if (!is.finite(sse)) return(list(log_target = -Inf, sse = Inf))
    list(log_target = -(alpha + n / 2) * log(beta + sse / 2), sse = sse)
  }
  current <- evaluate(state)
  if (!is.finite(current$log_target)) stop("Nonfinite initial conditional target.")
  n_save <- floor((n_iter - burn) / thin)
  draws <- matrix(NA_real_, n_save, 4L, dimnames = list(NULL, c("d", "phi", "theta", "sig2")))
  accepted <- setNames(integer(length(active)), active)
  saved <- 0L
  for (iteration in seq_len(n_iter)) {
    for (parameter in active) {
      candidate <- state
      candidate[parameter] <- state[parameter] + rnorm(1L, sd = proposal_sd[parameter])
      # Ordinary symmetric random walk with rejection outside support; no
      # truncation normalization or proposal-density correction is required.
      if (abs(candidate[parameter]) < bounds[parameter]) {
        proposed <- evaluate(candidate)
        if (is.finite(proposed$log_target) &&
            log(runif(1L)) < proposed$log_target - current$log_target) {
          state <- candidate
          current <- proposed
          accepted[parameter] <- accepted[parameter] + 1L
        }
      }
    }
    if (iteration > burn && (iteration - burn) %% thin == 0L) {
      saved <- saved + 1L
      sig2 <- 1 / rgamma(1L, shape = alpha + n / 2, rate = beta + current$sse / 2)
      draws[saved, ] <- c(state, sig2)
    }
    if (isTRUE(verbose) && iteration %% 1000L == 0L)
      message("Conditional ARFIMA iteration ", iteration, "/", n_iter)
  }
  list(d = draws[, "d"], phi = if (p == 1) draws[, "phi"] else numeric(0),
       theta = if (q == 1) draws[, "theta"] else numeric(0), sig2 = draws[, "sig2"],
       draws = draws, accept = accepted / n_iter, Y = Y,
       method = "conditional Gaussian ARFIMA with zero presample values (approximate)",
       p = p, q = q, prior = list(d = c(-0.5, 0.5), phi = c(-1, 1),
                                theta = c(-1, 1), sigma2_shape = alpha, sigma2_scale = beta),
       control = list(n_iter = n_iter, burn = burn, thin = thin, seed = seed,
                      d_sd = d_sd, arma_sd = arma_sd),
       warning = paste("Not the exact Gabor-Grazian filtered benchmark.",
                       "No unknown mean is fitted; pass mean-zero data or document centering.",
                       "Inspect convergence and initialization sensitivity before scientific use."))
}
