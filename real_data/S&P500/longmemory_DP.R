# Categorical sampler 
rcat <- function(prob) {
  (1:length(prob)) %*% rmultinom(1, 1, ifelse(prob > 0, prob, 0))
}

# Compute pi_h from V_h via stick-breaking
stick_break <- function(V) {
  H <- length(V)
  pi <- numeric(H)
  prod_term <- 1
  for (h in 1:H) {
    pi[h] <- V[h] * prod_term
    prod_term <- prod_term * (1 - V[h])
  }
  pi
}

longmemory_dp_sampler <- function(
    Y,                # response variable: log-periodogram (or smoothed log-periodogram) values Y_{j,n}
    omega,            # corresponding frequencies (on the trimmed grid)
    H       = 20,     # truncation level for DP (number of mixture components)
    runs    = 30000,  # total MCMC iterations
    burn    = 6000,   # burn-in
    thin    = 5,      # thinning
    a_V     = 1,      # Beta(a_V, b_V) hyperparameters for V_h
    b_V     = 2,
    a_sigma = 0.1,    # IG(a_sigma, b_sigma) prior for σ^2
    b_sigma = 0.1,
    a_xi    = 0.1,    # IG(a_xi, b_xi) prior for σ_xi^2
    b_xi    = 0.1,
    d_init  = 0,      # initial value for d
    tau_d   = 0.05,   # proposal sd for random-walk Metropolis on d
    display = 1000    # how often to print iteration counter
) {
  n <- length(Y)
  if (length(omega) != n) stop("Y and omega must have same length")
  
  ## --- Model: log f_d(ω) = log c_d - 2 d log |ω|
  ## For simplicity we set c_d = 1 here; it can be folded into ξ_h if desired.
  log_f_d <- function(d, omega) {
    -2 * d * log(abs(omega))
  }
  
  ## --- Storage sizes
  n_save <- floor((runs - burn) / thin)
  d_save     <- numeric(n_save)
  sigma2_save <- numeric(n_save)
  xi_save     <- matrix(NA, nrow = n_save, ncol = H)
  V_save      <- matrix(NA, nrow = n_save, ncol = H)
  Z_save      <- matrix(NA, nrow = n_save, ncol = n)  # optional: can be large and not saved if the series is long
  
  ## --- Initial values
  set.seed(123)
  d      <- d_init
  sigma2 <- 1
  sigma2_xi <- 1
  xi     <- rnorm(H, 0, sqrt(sigma2_xi))
  V      <- rbeta(H, a_V, b_V)
  pi_h   <- stick_break(V)
  Z      <- sample(1:H, n, replace = TRUE, prob = pi_h)
  
  ## --- Log-likelihood function in d (for MH update)
  loglik_d <- function(d, Y, omega, Z, xi, sigma2) {
    mu <- log_f_d(d, omega) + xi[Z]
    sum(dnorm(Y, mean = mu, sd = sqrt(sigma2), log = TRUE))
  }
  
  ## --- MCMC
  save_idx <- 0
  
  for (iter in 1:runs) {
    if (iter %% display == 0) cat("Iteration:", iter, "\n")
    
    ## 1. Update allocations Z_j
    log_f <- log_f_d(d, omega)
    for (j in 1:n) {
      log_prob <- numeric(H)
      for (h in 1:H) {
        log_prob[h] <- log(pi_h[h]) +
          dnorm(Y[j], mean = log_f[j] + xi[h], sd = sqrt(sigma2), log = TRUE)
      }
      # stabilise
      mmax <- max(log_prob)
      prob <- exp(log_prob - mmax)
      prob <- prob / sum(prob)
      Z[j] <- rcat(prob)
    }
    
    ## Component counts
    n_h <- tabulate(Z, nbins = H)
    
    ## 2. Update V_h and pi_h (stick-breaking weights)
    for (h in 1:H) {
      n_right <- if (h < H) sum(n_h[(h + 1):H]) else 0
      V[h] <- rbeta(1, a_V + n_h[h], b_V + n_right)
    }
    pi_h <- stick_break(V)
    
    ## 3. Update xi_h (atoms)
    log_f <- log_f_d(d, omega)
    for (h in 1:H) {
      idx <- which(Z == h)
      n_hh <- length(idx)
      if (n_hh > 0) {
        # posterior variance
        v_h <- 1 / (n_hh / sigma2 + 1 / sigma2_xi)
        # posterior mean
        m_h <- v_h * (1 / sigma2) * sum(Y[idx] - log_f[idx])
        xi[h] <- rnorm(1, mean = m_h, sd = sqrt(v_h))
      } else {
        # no data allocated: revert to prior
        xi[h] <- rnorm(1, mean = 0, sd = sqrt(sigma2_xi))
      }
    }
    
    ## 4. Update sigma_xi^2 (variance of atoms)
    # prior sigma_xi^2 ~ IG(a_xi, b_xi), xi_h ~ N(0, sigma_xi^2)
    shape_xi <- a_xi + H / 2
    rate_xi  <- b_xi + sum(xi^2) / 2
    sigma2_xi <- 1 / rgamma(1, shape = shape_xi, rate = rate_xi)
    
    ## 5. Update sigma^2 (noise variance)
    r <- Y - log_f_d(d, omega) - xi[Z]
    shape_sig <- a_sigma + n / 2
    rate_sig  <- b_sigma + sum(r^2) / 2
    sigma2 <- 1 / rgamma(1, shape = shape_sig, rate = rate_sig)
    
    ## 6. Update d via random-walk Metropolis
    d_prop <- d + tau_d * rnorm(1)
    if (d_prop > -1 && d_prop < 0.5) {
      ll_cur  <- loglik_d(d,      Y, omega, Z, xi, sigma2)
      ll_prop <- loglik_d(d_prop, Y, omega, Z, xi, sigma2)
      acc <- exp(ll_prop - ll_cur)  # prior is uniform on (-1,1/2)
      if (runif(1) < acc) d <- d_prop
    }
    
    ## --- Save draws
    if (iter > burn && ((iter - burn) %% thin == 0)) {
      save_idx <- save_idx + 1
      d_save[save_idx]       <- d
      sigma2_save[save_idx]  <- sigma2
      xi_save[save_idx, ]    <- xi
      V_save[save_idx, ]     <- V
      Z_save[save_idx, ]     <- Z
    }
  }
  
  list(
    d        = d_save,
    sigma2   = sigma2_save,
    xi       = xi_save,
    V        = V_save,
    Z        = Z_save,
    omega    = omega,
    Y        = Y,
    H        = H,
    a_V      = a_V,
    b_V      = b_V,
    a_sigma  = a_sigma,
    b_sigma  = b_sigma,
    a_xi     = a_xi,
    b_xi     = b_xi
  )
}

