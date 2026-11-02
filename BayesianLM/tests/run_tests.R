# Run from any directory with Rscript /path/to/project/tests/run_tests.R.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
root <- if (length(file_arg)) {
  dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
} else normalizePath(".")
source(file.path(root, "R", "longmemory_ksb.R"))
source(file.path(root, "R", "periodogram.R"))
source(file.path(root, "R", "simulate_arfima.R"))
source(file.path(root, "R", "arfima_benchmark.R"))

checks <- 0L
check <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
  checks <<- checks + 1L
  cat("PASS:", label, "\n")
}
near <- function(x, y, tolerance = 1e-10) {
  isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))
}
errors <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")

# Density identities, independent of the sampler implementation.
t <- c(0, .12, .5, 1)
v <- c(.4, .7)
knots <- c(.2, .8)
lambda <- 3
w <- ksb_weights(t, v, knots, lambda)
z1 <- v[1] * exp(-lambda * (t - knots[1])^2)
z2 <- v[2] * exp(-lambda * (t - knots[2])^2)
expected_w <- cbind(z1, (1 - z1) * z2, (1 - z1) * (1 - z2))
check(near(w, expected_w), "kernel stick weights and terminal residual")
check(near(rowSums(w), rep(1, length(t))), "weights sum to one at every covariate")
check(all(w >= 0), "weights nonnegative")
check(!near(w[1, ], w[4, ]), "weights vary across covariates")
check(near(as.numeric(ksb_weights(t, numeric(0), numeric(0), lambda)), rep(1, length(t))),
      "single-component truncation has unit weight")

p <- .27; mu <- 1.3; s1 <- .7; s2 <- 1.6
atom <- function(u) ksb_atom_density(u, p, mu, s1, s2)
check(abs(integrate(atom, -Inf, Inf)$value - 1) < 1e-7,
      "asymmetric atom integrates to one")
check(abs(integrate(function(u) u * atom(u), -Inf, Inf)$value) < 1e-7,
      "asymmetric atom has zero first moment")
u <- c(-3, -.2, 1, 5)
ps <- c(.2, .6, .4); mus <- c(1, -.7, .3)
s1s <- c(.5, 1.2, .8); s2s <- c(1.1, .6, 1.4)
direct <- vapply(seq_along(u), function(j) {
  sum(expected_w[j, ] * (ps * dnorm(u[j], mus, s1s) +
      (1 - ps) * dnorm(u[j], -ps * mus / (1 - ps), s2s)))
}, numeric(1))
ld <- ksb_log_density(u, t, v, knots, lambda, ps, mus, s1s, s2s)
check(near(exp(ld), direct), "marginal density equals explicit Gaussian sum")
check(all(is.finite(ksb_atom_density(c(-1e3, 1e3), p, mu, s1, s2, log = TRUE))),
      "log density remains finite in extreme tails")
pvalues <- c(.1, .3, .6, .9)
pi_target <- vapply(pvalues, function(p) ksb_pi_logtarget(p, numeric(), integer(), 0, 1, 1), numeric(1))
expected_pi <- dbeta(pvalues, 2, 2, log = TRUE) + log(pvalues) + log1p(-pvalues)
check(near(pi_target - pi_target[1], expected_pi - expected_pi[1]),
      "pi logit target uses Beta(2,2) with exactly one Jacobian")

# Raw FFT normalization and tapered Fourier-bin conventions.
n <- 128L; tt <- seq_len(n); k <- 7L
series <- cos(2 * pi * k * tt / n) + .31 * sin(2 * pi * 11 * tt / n)
raw <- periodogram_data(series, cutoff = 20, taper = "none")
check(near(raw$omega, 2 * pi * seq_len(20) / n), "raw Fourier frequencies")
check(near(raw$I[k], n / (8 * pi)), "raw periodogram normalization for unit cosine")
check(near(raw$x, log(4 * sin(raw$omega / 2)^2)), "exact log-periodogram regressor")
check(near(raw$t, raw$omega / max(raw$omega)), "conditional-density covariate normalization")
tapered <- periodogram_data(series, cutoff = 21, taper = "order1")
check(identical(as.integer(tapered$indices), seq.int(2L, 20L, by = 2L)),
      "order-one taper retains every second bin")
h <- 1 - exp(2i * pi * tt / n)
direct_I <- Mod(fft(h * (series - mean(series)), inverse = TRUE))^2 / (2 * pi * sum(Mod(h)^2))
check(near(tapered$I, direct_I[tapered$indices + 1L]), "tapered FFT and energy normalization")

# Construct a real sequence with exactly specified raw Fourier magnitudes.
# Its log-periodogram slope is exactly -d; this catches sign and double-log bugs.
target_d <- .3
freq <- 2 * pi * seq_len(n / 2 - 1L) / n
x <- log(4 * sin(freq / 2)^2)
fcoef <- complex(length = n)
mag <- sqrt(2 * pi * n * exp(-target_d * x))
fcoef[2:(n / 2)] <- mag
fcoef[(n / 2 + 2):n] <- rev(mag)
exact_series <- Re(fft(fcoef, inverse = TRUE) / n)
gph <- gph_estimate(exact_series, cutoff = 40, first = 5)
check(abs(gph$estimate - target_d) < 1e-9, "GPH memory sign and no double logarithm")
check(gph$se > 0 && gph$lower < gph$estimate && gph$upper > gph$estimate,
      "GPH reports nonzero asymptotic uncertainty")

# Small seeded integration run: computational smoke test, NOT convergence evidence.
omega <- 2 * pi * (1:24) / 128
Y <- -.8 - .2 * log(4 * sin(omega / 2)^2) + sin(1:24)
small_fit <- function(seed = 123L, ...) longmemory_ksb_sampler(
  Y, omega, H = 3, runs = 30, burn = 12, thin = 3, seed = seed, display = 0, ...)
a <- small_fit(); b <- small_fit()
check(length(a$d) == 6L && length(a$c) == 6L, "retained draw count")
check(all(is.finite(a$d)) && all(is.finite(a$c)), "finite regression draws")
check(all(a$d > -1 & a$d < .5) && all(abs(a$c) <= 50), "regression prior support")
check(near(a$d, b$d) && near(a$c, b$c), "same explicit seed reproduces draws")
cc <- small_fit(seed = 124L)
check(!near(a$d, cc$d) || !near(a$c, cc$c), "different seeds give different draws")
set.seed(234L); unseeded1 <- small_fit(seed = NULL); unseeded2 <- small_fit(seed = NULL)
check(!near(unseeded1$d, unseeded2$d) || !near(unseeded1$c, unseeded2$c),
      "seed=NULL does not reset the random number stream")
check(all(a$alpha > 0) && all(a$lambda > 0), "positive hyperparameter draws")
check(all(a$sigma1_sq >= .1 & a$sigma1_sq <= 50) &&
      all(a$sigma2_sq >= .1 & a$sigma2_sq <= 50), "variance bounds apply to variances")
check(all(abs(a$mu) <= 10) && all(a$pi > 0 & a$pi < 1), "atom prior support")
saved_ll <- sum(ksb_log_density(Y - a$c[1] + a$d[1] * a$x, a$t,
  a$V[1, ], a$Gamma[1, ], a$lambda[1], a$pi[1, ], a$mu[1, ],
  sqrt(a$sigma1_sq[1, ]), sqrt(a$sigma2_sq[1, ])))
check(near(saved_ll, a$log_pseudolik[1]), "saved pseudo-likelihood matches saved model state")
single <- longmemory_ksb_sampler(Y, omega, H = 1, runs = 6, burn = 2,
                               thin = 2, seed = 57, display = 0)
check(length(single$d) == 2 && ncol(single$V) == 0 && all(is.finite(single$d)),
      "single-atom sampler handles empty stick vectors")
check(errors(longmemory_ksb_sampler(Y[-1], omega, runs = 3, burn = 1)),
      "reject mismatched observations and frequencies")
check(errors(longmemory_ksb_sampler(Y, replace(omega, 1, 0), runs = 3, burn = 1)),
      "reject zero frequency")
check(errors(longmemory_ksb_sampler(Y, omega, runs = 3, burn = 3)), "reject invalid burn-in")

sim1 <- simulate_arfima(128, .2, burn = 512, seed = 19)
sim2 <- simulate_arfima(128, .2, burn = 512, seed = 19)
check(length(sim1) == 128 && all(is.finite(sim1)), "ARFIMA simulation length and finiteness")
check(near(sim1, sim2), "ARFIMA simulation seed reproducibility")
simt <- simulate_arfima(128, .2, innovations = "student", df = 5, burn = 512, seed = 19)
check(length(simt) == 128 && all(is.finite(simt)), "Student-t simulation smoke test")
white <- simulate_arfima(128, 0, burn = 0, seed = 38)
set.seed(38); expected_white <- rnorm(128)
check(near(white, expected_white), "ARFIMA(0,0,0) recovers innovations exactly")
arma <- simulate_arfima(128, 0, phi = .4, theta = .3, burn = 0, seed = 38)
expected_arma <- numeric(128)
for (j in seq_len(128)) {
  expected_arma[j] <- expected_white[j] + if (j > 1L) {
    .3 * expected_white[j - 1L] + .4 * expected_arma[j - 1L]
  } else 0
}
check(near(arma, expected_arma), "ARMA recursion and MA sign agree with manuscript")

check(near(arfima_fractional_difference(expected_white, 0), expected_white),
      "fractional difference at zero is identity")
check(near(arfima_fractional_difference(expected_white, 1),
           c(expected_white[1], diff(expected_white))),
      "fractional difference at one is first difference")
check(near(arfima_conditional_innovations(expected_arma, 0, .4, .3), expected_white),
      "conditional ARMA filter inverts independent recursion")
benchmark <- run_filtered_mcmc_arfima(expected_arma, p = 1, q = 1,
  n_iter = 24, burn = 12, thin = 3, seed = 89)
check(nrow(benchmark$draws) == 4 && all(is.finite(benchmark$draws)),
      "conditional benchmark finite saved draws")
check(all(abs(benchmark$d) < .5) && all(abs(benchmark$phi) < 1) &&
      all(abs(benchmark$theta) < 1) && all(benchmark$sig2 > 0),
      "conditional benchmark draws respect stated prior support")

cat(sprintf("\nAll %d checks passed. Short chains do not establish convergence.\n", checks))
