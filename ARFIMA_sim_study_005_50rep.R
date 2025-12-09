library(arfima)
source("longmemory_DP.R")   # contains longmemory_dp_sampler()
sourve("mcmc_simulation.R") # contains run_filtered_mcmc_arfima()

N <- 10000
lambda <- 2 * 0:N * pi / N
lambda <- lambda[1:(N/2)]
n <- floor(N/2)
lambdaj <- 2*pi/N * 1:n

m <- floor(1 + N^0.8)
l <- 1

# For diagnostics
y_mat    <- matrix(NA, ncol = length(l:m), nrow = 50)

H <- 20       
burn <- 6000
runs <- 30000
thin <- 5     

# Storing mixture components
y_mat  <- matrix(NA, ncol = length(l:m), nrow = 50)
xi_mat <- matrix(NA, ncol = H, row = 50)
V_mat  <- matrix(NA, ncol = H, row = 50)

# Storing the estimates of d
d_semi_vec   <- numeric(50)
d_gph_vec <- numeric(50)
d_filtered <- numeric(50)

# TRUE parameters for the ARFIMA simulation (adjust as needed)
d_true     <- 0.05
phi_true   <- 0
theta_true <- 0

# Optional: store allocations
z_mat <- matrix(NA, ncol = length(l:m), nrow = 50)

for(rep in 1:50)
{
  cat("Replication:", rep, "\n")
  set.seed(6533 + 10*rep)
  
  sim <- arfima.sim(N, model = list(phi = phi_true, theta = theta_true, dfrac = d_true))
  
  # Periodogram
  I  <- abs(fft(sim)/sqrt(N))^2
  P  <- (4/N) * I[1:(N/2)]
  logIj <- log(P[-1])            # remove zero frequency
  
  # trimmed grid
  Yj    <- logIj[l:m]
  omega <- lambdaj[l:m]
  
  y      <- Yj
  n_trim <- length(y)
  y_mat[rep, ] <- y
  
  # STEP 1: GPH estimator 
  d_gph_vec[rep] <- d_GPH(y, omega)

  # STEP 2: Bayesian semiparametric approach
  obj <- longmemory_dp_sampler(
    Y       = y,
    omega   = omega,
    H       = H,
    runs    = runs,
    burn    = burn,
    thin    = thin,
    a_V     = 1,
    b_V     = 2,
    a_sigma = 0.1,
    b_sigma = 0.1,
    a_xi    = 0.1,
    b_xi    = 0.1,
    d_init  = 0,
    tau_d   = 0.05,
    display = 1000
  )  
  
  # posterior mean of d under semiparametric model
  d_semi_vec[rep] <- mean(obj$d)
  
  # store mixture components (optional)
  xi_mat[rep, ] <- colMeans(obj$xi)
  V_mat[rep, ]  <- colMeans(obj$V)
  
  ## STEP 3: Filtered MCMC (Gabor & Grazian (2026))
  res_filt <- run_filtered_mcmc_arfima(
    Y       = sim,
    p       = 0,          #to change for other ARFIMA models
    q       = 0,         #to change for other ARFIMA models
    n_iter  = 30000,
    burn    = 6000,
    d_sd    = 0.05,
    verbose = TRUE
  )
  
  d_filtered[rep] = mean(res_filt$d)
  
  # ---- Save main posterior quantities ----
  
  # posterior mean of d
  d_vec[rep] <- mean(obj$d)
  
  # save Y for inspection
  y_mat[rep,] <- y
  
  # mixture components
  xi_mat[rep, ] <- colMeans(obj$xi)
  V_mat[rep, ]  <- colMeans(obj$V)
  
  # modal allocations Z_j
  for(j in 1:n_trim) {
    z_post <- table(obj$Z[, j])
    z_mat[rep, j] <- as.numeric(names(z_post)[which.max(z_post)])
  }
  
  save.image("BayesARFIMA-005.RData")
}

save.image("BayesARFIMA-005.RData")
