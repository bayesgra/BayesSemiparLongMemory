# ------------------------------------------------------------------
# Filtered MCMC for ARFIMA(p,d,q) as in Cohen-Gabor & Grazian (Stat&Comp 2025)
# ------------------------------------------------------------------
# Requires packages:
#   arfima    (for mle.arfima)
#   nsarfima  (for arfima.sim and tacvfARFIMA)
#   sarima    (for xarmaFilter)
#   extraDistr (for rinvgamma)
#   truncnorm  (for rtruncnorm, dtruncnorm)
#   tmvtnorm   (for rtmvnorm, dtmvnorm)
#
# Y:       observed time series
# p,q:     AR and MA orders (0 or 1 in this implementation)
# n_iter:  total MCMC iterations
# burn:    burn-in
# d_sd:    RW sd for d (truncated normal)
# alpha,beta: IG prior hyperparameters for sigma^2
# arma_sigm_reduction_factor: scale factor for ARMA proposal covariance
# verbose: print progress?
#
# Returns:
#   list(
#     d      = posterior draws of d (after burn),
#     phi    = posterior draws of phi (length 0 if p = 0),
#     theta  = posterior draws of theta (length 0 if q = 0),
#     sig2   = posterior draws of sigma^2,
#     accept = list(d = acc_rate_d, arma = acc_rate_arma),
#     Y      = original data
#   )
# ------------------------------------------------------------------

run_filtered_mcmc_arfima <- function(
    Y,
    p = 1, q = 1,
    n_iter = 5000,
    burn   = 1000,
    d_sd   = 0.025,
    alpha  = 1,
    beta   = 1,
    arma_sigm_reduction_factor = 1,
    verbose = TRUE
) {
  # ---- Checks ----
  if (p < 0 || q < 0 || p > 1 || q > 1)
    stop("This implementation currently supports p,q in {0,1} only.")
  
  n <- length(Y)
  
  # Helper: ARMA acf and ARFIMA(0,d,0) acf
  post_d_arfima0d0 <- function(Z, d, alpha, beta) {
    nZ <- length(Z)
    # ARFIMA(0,d,0) autocov via tacvfARFIMA (phi=0, theta=0, d=d)
    acf_vec <- nsarfima::tacvfARFIMA(phi = 0, theta = 0, dfrac = d, maxlag = nZ - 1)
    cov_mat <- toeplitz(acf_vec)
    Sig_inv <- solve(cov_mat)
    log_det <- determinant(cov_mat, logarithm = TRUE)$modulus
    
    Ev <- as.numeric(t(Z) %*% Sig_inv %*% Z)
    if (Ev < 0) return(-Inf)
    
    -0.5 * log_det - (nZ/2 + alpha) * log(Ev + 2 * beta)
  }
  
  post_arma <- function(U, phi, theta, alpha, beta) {
    nU <- length(U)
    # ARMA(p,q) autocovariance
    acf_vec <- ARMAacf(ar = if (p == 1) phi else numeric(0),
                       ma = if (q == 1) theta else numeric(0),
                       lag.max = nU - 1)
    cov_mat <- toeplitz(acf_vec)
    cov_mat[abs(cov_mat) < 1e-8] <- 0
    Sig_inv <- solve(cov_mat)
    log_det <- determinant(cov_mat, logarithm = TRUE)$modulus
    
    Ev <- as.numeric(t(U) %*% Sig_inv %*% U)
    if (Ev < 0) return(-Inf)
    
    -0.5 * log_det - (nU/2 + alpha) * log(Ev + 2 * beta)
  }
  
  # Full conditional for sigma^2 under ARFIMA(p,d,q)
  update_sig2 <- function(Y, d, phi, theta, alpha, beta) {
    acf_vec <- nsarfima::tacvfARFIMA(phi = if (p == 1) phi else 0,
                                     theta = if (q == 1) -theta else 0,
                                     dfrac = d,
                                     maxlag = n - 1)
    cov_mat <- toeplitz(acf_vec)
    Sig_inv <- solve(cov_mat)
    Ev <- as.numeric(t(Y) %*% Sig_inv %*% Y)
    shape <- n/2 + alpha
    rate  <- (Ev + 2 * beta) / 2
    extraDistr::rinvgamma(1, shape = shape, rate = rate)
  }
  
  # ---- Initial values from ARFIMA MLE ----
  # Use arfima::mle.arfima for starting values
  mle_fit <- arfima::mle.arfima(Y, p = p, q = q, d.range = c(0, 0.5), incl.mean = FALSE)
  par0    <- mle_fit$pars
  
  d0     <- as.numeric(par0["d"])
  sig20  <- as.numeric(par0["sig2"])
  phi0   <- if (p == 1) {
    # try to pick the AR coefficient (often called "ar1" or similar)
    as.numeric(par0[grep("^ar", names(par0))[1]])
  } else 0
  theta0 <- if (q == 1) {
    as.numeric(par0[grep("^ma", names(par0))[1]])
  } else 0
  
  # ARMA MLE on differenced series to get proposal covariance
  U0 <- fracdiff::diffseries(Y, d = d0)[-1]
  if (p == 0 && q == 0) {
    arma_mle  <- numeric(0)
    arma_cov  <- matrix(0, 0, 0)
  } else {
    arma_fit  <- stats::arima(U0, order = c(p, 0, q))
    arma_coef <- arma_fit$coef
    arma_names <- character(0)
    if (p == 1) arma_names <- c(arma_names, "ar1")
    if (q == 1) arma_names <- c(arma_names, "ma1")
    arma_mle <- arma_coef[arma_names]
    names(arma_mle) <- c(if (p == 1) "phi", if (q == 1) "theta")
    
    arma_cov <- arma_fit$var.coef[arma_names, arma_names, drop = FALSE] /
      arma_sigm_reduction_factor
  }
  
  # ---- Storage ----
  d_chain     <- numeric(n_iter)
  phi_chain   <- numeric(n_iter)
  theta_chain <- numeric(n_iter)
  sig2_chain  <- numeric(n_iter)
  
  # initial values
  d_curr     <- d0
  phi_curr   <- phi0
  theta_curr <- theta0
  sig2_curr  <- sig20
  
  # initial posteriors
  # whiten Y with current ARMA for d-update
  Z_curr <- if (p == 0 && q == 0) Y else {
    sarima::xarmaFilter(x = Y, model = list(ar = phi_curr, ma = theta_curr), whiten = TRUE)
  }
  post_d_curr <- post_d_arfima0d0(Z_curr, d_curr, alpha, beta)
  
  # ARMA posterior given differenced series
  U_curr <- fracdiff::diffseries(Y, d = d_curr)[-1]
  post_arma_curr <- if (p == 0 && q == 0) 0 else {
    post_arma(U_curr, phi_curr, theta_curr, alpha, beta)
  }
  
  accept_d    <- 0
  accept_arma <- 0
  
  # ---- MCMC loop ----
  for (iter in 1:n_iter) {
    if (verbose && iter %% 500 == 0) {
      cat("Iteration", iter, "d =", round(d_curr,3),
          "phi =", round(phi_curr,3),
          "theta =", round(theta_curr,3), "\n")
    }
    
    ## -------------------------
    ## 1. Update d (filtered part)
    ## -------------------------
    if (p == 0 && q == 0) {
      # If no ARMA, Z = Y (no whitening)
      Z_curr <- Y
    } else {
      Z_curr <- sarima::xarmaFilter(
        x     = Y,
        model = list(ar = phi_curr, ma = theta_curr),
        whiten = TRUE
      )
    }
    
    d_prop <- truncnorm::rtruncnorm(1, mean = d_curr, sd = d_sd,
                                    a = -0.5, b = 0.5)
    post_d_prop <- post_d_arfima0d0(Z_curr, d_prop, alpha, beta)
    
    # proposal densities for truncated normal
    q_curr_given_prop <- truncnorm::dtruncnorm(
      x = d_curr, a = -0.5, b = 0.5,
      mean = d_prop, sd = d_sd, log = TRUE
    )
    q_prop_given_curr <- truncnorm::dtruncnorm(
      x = d_prop, a = -0.5, b = 0.5,
      mean = d_curr, sd = d_sd, log = TRUE
    )
    
    log_acc_d <- post_d_prop + q_curr_given_prop -
      post_d_curr - q_prop_given_curr
    if (is.finite(log_acc_d) && log(runif(1)) < log_acc_d) {
      d_curr      <- d_prop
      post_d_curr <- post_d_prop
      accept_d    <- accept_d + 1
    }
    
    ## -------------------------
    ## 2. Update ARMA parameters (phi, theta)
    ## -------------------------
    if (p > 0 || q > 0) {
      # Differenced series for ARMA
      U_curr <- fracdiff::diffseries(Y, d = d_curr)[-1]
      k <- p + q
      arma_now <- c(if (p == 1) phi_curr, if (q == 1) theta_curr)
      names(arma_now) <- c(if (p == 1) "phi", if (q == 1) "theta")
      
      arma_prop <- as.numeric(
        tmvtnorm::rtmvnorm(
          1,
          mean  = arma_now,
          sigma = arma_cov,
          lower = rep(-1, k),
          upper = rep( 1, k)
        )
      )
      names(arma_prop) <- names(arma_now)
      
      phi_prop   <- if (p == 1) arma_prop["phi"]   else 0
      theta_prop <- if (q == 1) arma_prop["theta"] else 0
      
      post_arma_prop <- post_arma(U_curr, phi_prop, theta_prop, alpha, beta)
      
      # proposal densities (truncated MVN)
      log_q_curr_given_prop <- tmvtnorm::dtmvnorm(
        x     = arma_now,
        mean  = arma_prop,
        sigma = arma_cov,
        lower = rep(-1, k),
        upper = rep( 1, k),
        log   = TRUE
      )
      log_q_prop_given_curr <- tmvtnorm::dtmvnorm(
        x     = arma_prop,
        mean  = arma_now,
        sigma = arma_cov,
        lower = rep(-1, k),
        upper = rep( 1, k),
        log   = TRUE
      )
      
      log_acc_arma <- post_arma_prop + log_q_curr_given_prop -
        post_arma_curr - log_q_prop_given_curr
      
      if (is.finite(log_acc_arma) && log(runif(1)) < log_acc_arma) {
        phi_curr        <- phi_prop
        theta_curr      <- theta_prop
        post_arma_curr  <- post_arma_prop
        accept_arma     <- accept_arma + 1
      }
    }
    
    ## -------------------------
    ## 3. Update sigma^2
    ## -------------------------
    sig2_curr <- update_sig2(Y, d_curr, phi_curr, theta_curr, alpha, beta)
    
    ## -------------------------
    ## 4. Store
    ## -------------------------
    d_chain[iter]     <- d_curr
    phi_chain[iter]   <- phi_curr
    theta_chain[iter] <- theta_curr
    sig2_chain[iter]  <- sig2_curr
  }
  
  keep <- (burn+1):n_iter
  
  list(
    d      = d_chain[keep],
    phi    = if (p == 1) phi_chain[keep]   else numeric(0),
    theta  = if (q == 1) theta_chain[keep] else numeric(0),
    sig2   = sig2_chain[keep],
    accept = list(
      d    = accept_d    / n_iter,
      arma = accept_arma / n_iter
    ),
    Y      = Y
  )
}
