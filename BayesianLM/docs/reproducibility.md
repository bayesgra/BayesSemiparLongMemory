# Reproducibility and validation scope

## What this project can establish

The source implements the current manuscript's finite, mean-zero kernel
stick-breaking model and provides repeatable computational workflows. The tests
check mathematical and software identities and short chain behavior. A successful
test run is not evidence that a research-length chain has converged, that the
paper's empirical claims hold for this model, or that the full published study
has been reproduced.

## Explicit conventions where the manuscript is incomplete

1. **Kernel:** `K(t,Gamma,lambda) = exp(-lambda*(t-Gamma)^2)`. Here larger `lambda`
   gives a narrower kernel. The manuscript says squared exponential but does not
   specify its formula. A length-scale formula such as
   `exp(-(t-Gamma)^2/(2*lambda^2))`, with the same Gamma prior on `lambda`, gives a
   different prior and is not interchangeable with this choice.
2. **Order-one taper:** `h[s] = 1-exp(2*pi*i*s/n)`, with the positive-exponent
   Fourier sum. Normalize the squared modulus by `2*pi*sum(abs(h)^2)`. The raw
   periodogram is available separately. This explicit convention must be retained
   when comparing another implementation of the tapered transform.
   This is the convention in [Faÿ, equation (2.1)](https://arxiv.org/pdf/0807.5096).
3. **Low-frequency cutoff:** `M=floor(n^0.8)`. For the tapered analysis retain
   indices `2,4,...,2*floor(M/2)`; use their actual frequencies for the regression
   and normalize by the largest retained frequency for the mixture covariate.
4. **MCMC settings:** the implementation defaults to 30,000 iterations, 6,000
   burn-in iterations and thinning every 5 iterations. These are computational
   choices; the supplied manuscript does not report settings for the current
   sampler. Short demonstrations use much smaller counts.
5. **ARFIMA generation:** finite zero-start convolution with discarded warmup is
   an approximation to stationary sampling. The default warmup is
   `max(5000,5*n)`; check sensitivity especially near `d=0.5`. The MA sign is
   positive, matching `(1-phi*B)(1-B)^d X = (1+theta*B) epsilon`.
6. **Frequentist comparator:** the base-R routine uses raw log-periodogram
   regression starting at Fourier index 5 in simulations (index 1 in applications), with
   asymptotic log-error variance
   `pi^2/6`. It is explicitly labeled as such. It has not been certified identical
   to the paper's `LongMemoryTS::gph` procedure or its definition of `ell=5`.

## Full simulation design from Section 5

| Dimension | Values |
|---|---|
| Series length | 5,000; 10,000; 20,000 |
| Memory `d` | 0.05; 0.10; 0.20; 0.30; 0.40; 0.49 |
| AR/MA parameters | (0,0); (0.8,0.1); (0.5,0.1); (0.1,0.5); (0.1,0.8) |
| Replicates | 200 per setting |
| Primary innovations | Standard normal |
| Heavy-tail extension | Standardized Student-t with 5 degrees of freedom; ARFIMA(0,d,0) |
| Mixture truncation | 30 atoms |

There are 30 combinations per sample size and 90 combinations across the three
sample sizes. The initial cutoff counts are 910, 1,584 and 2,759; tapered retained
counts are 455, 792 and 1,379. The Gaussian study therefore requires 18,000 fits;
the heavy-tail extension adds 3,600 fits. Small demonstration runs are not this
full study.

## Before reporting results

Use independent seeds and overdispersed starting values for multiple chains.
Inspect traces and mixing of `c`, `d`, `lambda`, `alpha`, and the log
pseudo-likelihood; assess effective sample sizes and between-chain agreement.
Mixture labels may be poorly identified, so componentwise comparisons require
care. Inspect acceptance rates to guide proposal tuning; acceptance rates alone
are insufficient. Check sensitivity to truncation `H`, simulation warmup,
frequency cutoff and the stated kernel convention. Use chain diagnostics before
summarizing posterior quantiles or coverage.

The stored `terminal_weight_max` is the maximum terminal weight over **observed**
covariates for each saved draw, not the supremum over all of [0,1]. For a wider
check, evaluate `ksb_weights` on a dense covariate grid using saved parameters.

## Known replication gaps

- The supplied PDF references supplementary appendices, but those appendices are
  not included. Nile and temperature preprocessing cannot be recovered uniquely
  from the manuscript alone.
- The S&P application describes a full raw-frequency count, whereas the
  simulation section uses tapered and thinned ordinates. Application preprocessing
  must be stated explicitly with results.
- The optional `BayesConditional` comparator replaces the problematic cached
  target with a coherent zero-presample conditional likelihood. It does not
  implement the cited filtered MCMC method or reproduce the old `BayesParam`
  results; see `conditional-benchmark.md`.
- The old tables and figures were not recomputed for the new model. New outputs
  must be distinguished from those numerical claims.
- The original CSV data are bundled as supplied. Their inclusion does not grant
  redistribution rights or replace checking the data providers' attribution and
  license requirements before making a repository public.

## Software and provenance

The core model and workflows use base R; scripts should record `sessionInfo()`
alongside generated results. The author should choose a code license and verify
the bundled data's redistribution conditions before publishing. No license is
assigned on the author's behalf.
