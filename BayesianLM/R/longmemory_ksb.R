# Mean-zero kernel stick-breaking log-periodogram regression.
# Base R only. Normal priors below use VARIANCES 100 and 10.
# Squared-exponential convention: K(t,g,lambda)=exp(-lambda*(t-g)^2).
# The manuscript does not specify its bandwidth parameterization explicitly.

ksb_logadd <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

ksb_weights <- function(t, V, Gamma, lambda, log = FALSE) {
  if (!is.numeric(t) || !length(t) || any(!is.finite(t)) || any(t < 0 | t > 1) ||
      length(V) != length(Gamma) || any(!is.finite(V)) || any(V <= 0 | V >= 1) ||
      any(!is.finite(Gamma)) || any(Gamma < 0 | Gamma > 1) ||
      length(lambda) != 1L || !is.finite(lambda) || lambda <= 0)
    stop("Invalid kernel stick-breaking arguments")
  H <- length(V) + 1L
  lw <- matrix(0, length(t), H)
  remainder <- rep(0, length(t))
  if (H > 1L) for (h in seq_len(H - 1L)) {
    lz <- log(V[h]) - lambda * (t - Gamma[h])^2
    lw[, h] <- remainder + lz
    remainder <- remainder + log1p(-exp(lz))
  }
  lw[, H] <- remainder # all remaining mass, not a further random stick
  if (log) lw else exp(lw)
}

ksb_atom_density <- function(u, pi, mu, sigma1, sigma2, log = FALSE) {
  if (any(!is.finite(c(pi, mu, sigma1, sigma2))) || any(pi <= 0 | pi >= 1) ||
      any(sigma1 <= 0 | sigma2 <= 0)) stop("Invalid atom parameters")
  a <- log(pi) + dnorm(u, mu, sigma1, log = TRUE)
  b <- log1p(-pi) + dnorm(u, -pi * mu / (1 - pi), sigma2, log = TRUE)
  ans <- ksb_logadd(a, b)
  if (log) ans else exp(ans)
}

ksb_log_density <- function(u, t, V, Gamma, lambda, pi, mu, sigma1, sigma2) {
  H <- length(V) + 1L
  if (length(u) != length(t) || any(!is.finite(u)) ||
      any(lengths(list(pi, mu, sigma1, sigma2)) != H)) stop("Incompatible density dimensions")
  lw <- ksb_weights(t, V, Gamma, lambda, log = TRUE)
  ans <- rep(-Inf, length(u))
  for (h in seq_len(H)) ans <- ksb_logadd(ans, lw[, h] +
    ksb_atom_density(u, pi[h], mu[h], sigma1[h], sigma2[h], log = TRUE))
  ans
}

# Inverse-CDF truncated normal, switching tails to avoid cancellation.
ksb_rtnorm <- function(mean, sd, lower, upper) {
  a <- (lower - mean) / sd; b <- (upper - mean) / sd
  if (a > 0) {
    lo <- pnorm(b, lower.tail = FALSE); hi <- pnorm(a, lower.tail = FALSE)
    z <- qnorm(runif(1, lo, hi), lower.tail = FALSE)
  } else {
    lo <- pnorm(a); hi <- pnorm(b)
    z <- qnorm(runif(1, lo, hi))
  }
  if (!is.finite(z) || hi <= lo) stop("Truncated normal underflow; inspect scaling and chain state")
  pmin(upper, pmax(lower, mean + sd * z))
}

ksb_rig <- function(shape, rate, lower = 0.1, upper = 50) {
  # Inverse variance has a Gamma(shape, rate) distribution.
  a <- 1 / upper; b <- 1 / lower
  if (a > shape / rate) {
    lo <- pgamma(b, shape, rate = rate, lower.tail = FALSE)
    hi <- pgamma(a, shape, rate = rate, lower.tail = FALSE)
    if (hi <= lo || !is.finite(hi)) stop("Truncated inverse-gamma underflow")
    precision <- qgamma(runif(1, lo, hi), shape, rate = rate, lower.tail = FALSE)
  } else {
    lo <- pgamma(a, shape, rate = rate)
    hi <- pgamma(b, shape, rate = rate)
    if (hi <= lo || !is.finite(hi)) stop("Truncated inverse-gamma underflow")
    precision <- qgamma(runif(1, lo, hi), shape, rate = rate)
  }
  if (!is.finite(precision) || precision <= 0) stop("Invalid inverse-gamma draw")
  1 / precision
}

# Complete atom-conditional target on logit(pi), including its Jacobian.
ksb_pi_logtarget <- function(p, u, side, mu, sigma1_sq, sigma2_sq) {
  if (!is.finite(p) || p <= 0 || p >= 1) return(-Inf)
  2 * log(p) + 2 * log1p(-p) + sum(ifelse(side == 1L, log(p), log1p(-p))) +
    sum(dnorm(u, ifelse(side == 1L, mu, -p * mu / (1 - p)),
              sqrt(ifelse(side == 1L, sigma1_sq, sigma2_sq)), log = TRUE))
}

longmemory_ksb_sampler <- function(Y, omega, H = 30, runs = 30000,
    burn = 6000, thin = 5, seed = NULL, display = 0,
    d_init = NULL, c_init = NULL,
    proposal = list(pi = 0.35, V = 0.4, lambda = 0.2)) {
  is_int <- function(z, low) is.numeric(z) && length(z) == 1L &&
    is.finite(z) && z >= low && z == floor(z)
  if (!is.numeric(Y) || !is.numeric(omega) || length(Y) < 3L || length(Y) != length(omega) ||
      any(!is.finite(Y)) || any(!is.finite(omega)) || any(omega <= 0 | omega > base::pi) ||
      any(diff(omega) <= 0)) stop("Y and omega must be finite equal-length vectors; omega increasing in (0, pi]")
  if (!is_int(H, 1) || !is_int(runs, 1) || !is_int(burn, 0) || !is_int(thin, 1) ||
      !is_int(display, 0) || burn >= runs || floor((runs - burn) / thin) < 1L)
    stop("Invalid H, runs, burn, thin, or display")
  defaults <- list(pi = 0.35, V = 0.4, lambda = 0.2)
  if (!is.list(proposal) || is.null(names(proposal)) || any(!names(proposal) %in% names(defaults)))
    stop("proposal must be a named list with pi, V, and/or lambda")
  defaults[names(proposal)] <- proposal; proposal <- defaults
  if (any(!vapply(proposal, function(z) is.numeric(z) && length(z) == 1L &&
      is.finite(z) && z > 0, logical(1)))) stop("Proposal scales must be positive")
  if (!is.null(seed)) {
    if (!is_int(seed, 0)) stop("seed must be a nonnegative integer or NULL")
    set.seed(seed)
  }
  n <- length(Y); x <- log(4 * sin(omega / 2)^2); t <- omega / max(omega)
  X <- cbind(1, -x)
  initial <- lm.fit(X, Y)$coefficients
  d <- if (is.null(d_init)) max(-0.95, min(0.45, initial[2])) else d_init
  c0 <- if (is.null(c_init)) max(-49, min(49, mean(Y + d * x))) else c_init
  if (length(d) != 1L || !is.finite(d) || d <= -1 || d >= 0.5 ||
      length(c0) != 1L || !is.finite(c0) || abs(c0) > 50) stop("Initial regression parameters outside prior support")
  pi_atom <- rep(0.5, H); mu <- rnorm(H, 0, 0.3)
  s1 <- s2 <- rep(1, H); alpha <- 2; lambda <- 2
  V <- if (H > 1L) rbeta(H - 1L, 1, alpha) else numeric()
  Gamma <- if (H > 1L) runif(H - 1L) else numeric()
  z <- rep(1L, n); side <- rep(1L, n)
  ns <- floor((runs - burn) / thin)
  out <- list(c = numeric(ns), d = numeric(ns), alpha = numeric(ns), lambda = numeric(ns),
    pi = matrix(NA_real_, ns, H), mu = matrix(NA_real_, ns, H),
    sigma1_sq = matrix(NA_real_, ns, H), sigma2_sq = matrix(NA_real_, ns, H),
    V = matrix(NA_real_, ns, H - 1L), Gamma = matrix(NA_real_, ns, H - 1L),
    log_pseudolik = numeric(ns), terminal_weight_max = numeric(ns))
  accepted <- c(pi = 0, V = 0, Gamma = 0, lambda = 0, regression = 0)
  attempted <- accepted
  accept <- function(logratio) is.finite(logratio) && log(runif(1)) < logratio
  saved <- 0L
  for (iter in seq_len(runs)) {
    u <- Y - c0 + d * x
    lw <- ksb_weights(t, V, Gamma, lambda, log = TRUE)
    lp <- matrix(0, n, 2L * H)
    for (h in seq_len(H)) {
      lp[, h] <- lw[, h] + log(pi_atom[h]) + dnorm(u, mu[h], sqrt(s1[h]), log = TRUE)
      lp[, H + h] <- lw[, h] + log1p(-pi_atom[h]) +
        dnorm(u, -pi_atom[h] * mu[h] / (1 - pi_atom[h]), sqrt(s2[h]), log = TRUE)
    }
    # Gumbel-max sampling avoids rowwise normalization and probability underflow.
    choices <- max.col(lp - log(-log(matrix(runif(n * 2L * H), n))), ties.method = "first")
    z <- (choices - 1L) %% H + 1L; side <- ifelse(choices <= H, 1L, 2L)

    # Joint c,d Gaussian conditional, rejected outside the rectangular prior support.
    means <- ifelse(side == 1L, mu[z], -pi_atom[z] * mu[z] / (1 - pi_atom[z]))
    variances <- ifelse(side == 1L, s1[z], s2[z])
    precision <- crossprod(X, X / variances) + diag(c(1 / 100, 0))
    U <- chol(precision)
    center <- backsolve(U, forwardsolve(t(U), crossprod(X, (Y - means) / variances)))
    attempted["regression"] <- attempted["regression"] + 1
    for (attempt in seq_len(1000L)) {
      beta <- as.numeric(center + backsolve(U, rnorm(2)))
      if (abs(beta[1]) <= 50 && beta[2] > -1 && beta[2] < 0.5) {
        c0 <- beta[1]; d <- beta[2]
        accepted["regression"] <- accepted["regression"] + 1
        break
      }
    }
    u <- Y - c0 + d * x

    for (h in seq_len(H)) {
      idx <- which(z == h); ui <- u[idx]; ss <- side[idx]
      # Logit random walk: Beta(2,2) plus log-Jacobian gives 2 log p + 2 log(1-p).
      target_pi <- function(p) ksb_pi_logtarget(p, ui, ss, mu[h], s1[h], s2[h])
      pnew <- plogis(qlogis(pi_atom[h]) + rnorm(1, 0, proposal$pi))
      attempted["pi"] <- attempted["pi"] + 1
      if (accept(target_pi(pnew) - target_pi(pi_atom[h]))) {
        pi_atom[h] <- pnew; accepted["pi"] <- accepted["pi"] + 1
      }
      factors <- ifelse(ss == 1L, 1, -pi_atom[h] / (1 - pi_atom[h]))
      vv <- ifelse(ss == 1L, s1[h], s2[h])
      prec <- 1 / 10 + sum(factors^2 / vv)
      mu[h] <- ksb_rtnorm(sum(factors * ui / vv) / prec, sqrt(1 / prec), -10, 10)
      residual <- ui - factors * mu[h]
      for (k in 1:2) {
        rr <- residual[ss == k]
        variance <- ksb_rig(2 + length(rr) / 2, 2 + sum(rr^2) / 2)
        if (k == 1L) s1[h] <- variance else s2[h] <- variance
      }
    }

    if (H > 1L) for (h in seq_len(H - 1L)) {
      at <- z == h; right <- z > h
      stick_ll <- function(v, g) {
        lz <- log(v) - lambda * (t - g)^2
        sum(lz[at]) + sum(log1p(-exp(lz[right])))
      }
      vnew <- plogis(qlogis(V[h]) + rnorm(1, 0, proposal$V))
      attempted["V"] <- attempted["V"] + 1
      if (vnew > 0 && vnew < 1) {
        # Beta(1,alpha) density and logit Jacobian.
        ratio <- stick_ll(vnew, Gamma[h]) - stick_ll(V[h], Gamma[h]) +
          log(vnew) - log(V[h]) + alpha * (log1p(-vnew) - log1p(-V[h]))
        if (accept(ratio)) { V[h] <- vnew; accepted["V"] <- accepted["V"] + 1 }
      }
      gnew <- runif(1)
      attempted["Gamma"] <- attempted["Gamma"] + 1
      if (accept(stick_ll(V[h], gnew) - stick_ll(V[h], Gamma[h]))) {
        Gamma[h] <- gnew; accepted["Gamma"] <- accepted["Gamma"] + 1
      }
    }
    alpha <- rgamma(1, shape = 2 + H - 1L, rate = 1 - sum(log1p(-V)))
    target_lambda <- function(lam) {
      if (!is.finite(lam) || lam <= 0) return(-Inf)
      logw <- ksb_weights(t, V, Gamma, lam, log = TRUE)
      sum(logw[cbind(seq_len(n), z)]) + 2 * log(lam) - lam
    }
    lnew <- exp(log(lambda) + rnorm(1, 0, proposal$lambda))
    attempted["lambda"] <- attempted["lambda"] + 1
    if (accept(target_lambda(lnew) - target_lambda(lambda))) {
      lambda <- lnew; accepted["lambda"] <- accepted["lambda"] + 1
    }
    if (iter > burn && (iter - burn) %% thin == 0L) {
      saved <- saved + 1L
      out$c[saved] <- c0; out$d[saved] <- d
      out$alpha[saved] <- alpha; out$lambda[saved] <- lambda
      out$pi[saved, ] <- pi_atom; out$mu[saved, ] <- mu
      out$sigma1_sq[saved, ] <- s1; out$sigma2_sq[saved, ] <- s2
      out$V[saved, ] <- V; out$Gamma[saved, ] <- Gamma
      out$log_pseudolik[saved] <- sum(ksb_log_density(Y - c0 + d * x, t, V, Gamma,
        lambda, pi_atom, mu, sqrt(s1), sqrt(s2)))
      out$terminal_weight_max[saved] <- max(ksb_weights(t, V, Gamma, lambda)[, H])
    }
    if (display > 0L && iter %% display == 0L) message("Iteration ", iter, " / ", runs)
  }
  out$Y <- Y; out$omega <- omega; out$x <- x; out$t <- t; out$H <- H
  out$settings <- list(runs = runs, burn = burn, thin = thin, seed = seed, proposal = proposal,
    kernel = "exp(-lambda*(t-Gamma)^2)", d_bounds = c(-1, 0.5),
    c_variance = 100, c_bounds = c(-50, 50), mu_variance = 10, mu_bounds = c(-10, 10),
    variance_bounds = c(0.1, 50), pi_beta = c(2, 2), variance_ig = c(2, 2),
    alpha_gamma_shape_rate = c(2, 1), lambda_gamma_shape_rate = c(2, 1))
  rates <- accepted / attempted
  rates[attempted == 0] <- NA_real_
  out$diagnostics <- list(acceptance = rates, accepted = accepted, attempted = attempted,
    terminal_weight_mean = mean(out$terminal_weight_max),
    note = "Acceptance and finite traces do not establish convergence; compare independent chains and H values.")
  out$last_allocations <- list(component = z, normal_component = side)
  class(out) <- "longmemory_ksb_fit"
  out
}
