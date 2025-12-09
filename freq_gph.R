d_GPH <- function(I_per, omega, l = 1, m = floor(length(I_per)^0.8)) {
  # I_per = periodogram values at Fourier frequencies (excluding zero)
  # omega = corresponding frequencies
  # use only frequencies l,...,m
  
  idx <- l:m
  y   <- log(I_per[idx])
  
  x   <- log(4 * (sin(omega[idx] / 2))^2)
  
  fit <- lm(y ~ x)
  
  d_hat <- -coef(fit)[2]
  return(d_hat)
}

