# The order-one taper convention is h_s = 1 - exp(2*pi*i*s/n).
# Use a positive-exponent Fourier sum, matching Faÿ (2010), Eq. (2.1).
# The paper does not spell out the coefficient convention; see docs.
periodogram_data <- function(series, cutoff = floor(length(series)^0.8),
                             taper = c("order1", "none"), stride = NULL,
                             first = NULL, demean = TRUE) {
  taper <- match.arg(taper)
  if (!is.numeric(series) || length(series) < 8L || any(!is.finite(series)))
    stop("series must contain at least eight finite numeric observations")
  n <- length(series)
  if (is.null(stride)) stride <- if (taper == "order1") 2L else 1L
  if (is.null(first)) first <- if (taper == "order1") 2L else 1L
  integer_scalar <- function(z) length(z) == 1L && is.finite(z) && z == floor(z)
  if (!integer_scalar(cutoff) || cutoff < 1L || cutoff > floor(n / 2))
    stop("cutoff must be an integer between 1 and floor(n/2)")
  if (!integer_scalar(stride) || stride < 1L || !integer_scalar(first) ||
      first < 1L || first > cutoff) stop("invalid first frequency or stride")
  if (length(demean) != 1L || is.na(demean) || !is.logical(demean))
    stop("demean must be TRUE or FALSE")
  series <- as.numeric(series)
  if (demean) series <- series - mean(series)
  h <- if (taper == "order1") 1 - exp(2 * pi * 1i * seq_len(n) / n) else rep(1, n)
  indices <- seq.int(first, cutoff, by = stride)
  omega <- 2 * pi * indices / n
  I <- Mod(fft(series * h, inverse = TRUE)[indices + 1L])^2 / (2 * pi * sum(Mod(h)^2))
  if (any(!is.finite(I) | I <= 0))
    stop("retained periodogram includes nonpositive/nonfinite values; no flooring applied")
  list(Y = log(I), omega = omega, x = log(4 * sin(omega / 2)^2),
       t = omega / max(omega), indices = indices, I = I, n = n,
       cutoff = cutoff, taper = taper, stride = stride, first = first,
       normalization = "squared FFT modulus / (2*pi*sum(abs(taper)^2))")
}

# Classical raw-periodogram regression; model mean c - d*x.
# This does not call LongMemoryTS, and is explicitly labeled as a base-R comparator.
gph_estimate <- function(series, cutoff = floor(length(series)^0.8), first = 5L) {
  dat <- periodogram_data(series, cutoff, taper = "none", stride = 1L, first = first)
  if (length(dat$Y) < 3L) stop("GPH requires at least three retained ordinates")
  fit <- lm(dat$Y ~ dat$x)
  estimate <- -unname(coef(fit)[2L])
  se <- sqrt((pi^2 / 6) / sum((dat$x - mean(dat$x))^2))
  list(estimate = estimate, se = se, lower = estimate - qnorm(.975) * se,
       upper = estimate + qnorm(.975) * se, fit = fit, data = dat,
       method = "base-R raw-periodogram GPH, asymptotic variance pi^2/6")
}
