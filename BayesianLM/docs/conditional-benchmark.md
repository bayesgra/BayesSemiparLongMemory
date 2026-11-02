# Optional conditional ARFIMA benchmark

`R/arfima_benchmark.R` supplies a transparent base-R comparator. The legacy entry-point name `run_filtered_mcmc_arfima()` is retained for migration, but the method is explicitly **not** an implementation or numerical reproduction of the filtered MCMC method cited in the paper.

The model is `(1-phi B)(1-B)^d Y_t=(1+theta B)epsilon_t`, with all presample series and innovation values set to zero. This defines a conditional finite-past Gaussian likelihood and approximates the stationary ARFIMA likelihood. It is a single coherent target for all parameter updates, unlike the stale and incompatible cached conditional targets in the old script.

The fractional differencing coefficients are `w_0=1`, `w_k=w_{k-1}(k-1-d)/k`. The first n elements of the convolution are computed through zero-padded FFT multiplication. After AR filtering, the innovations solve `epsilon_t=z_t-theta*epsilon_{t-1}`. This implements the paper's plus-theta moving-average convention. The transformation from the observations to innovations is lower triangular with unit diagonal, so no parameter-dependent determinant is omitted.

Given sum of squared innovations S, integrating `sigma^2~IG(a,b)` gives a target proportional to `(b+S/2)^(-a-n/2)`. Defaults are `a=b=0.1`. Uniform priors are used on `d in (-0.5,0.5)` and the active AR/MA coefficients in `(-1,1)`, for known p,q in {0,1}. The narrower d range is specific to this conventional stationary/invertible comparator; the semiparametric method retains the manuscript's `(-1,0.5)` range. Each random-walk move evaluates the full current target under one common conditioning convention. Proposals outside support are rejected, so symmetric normal proposal terms cancel. Variance draws use their exact inverse-gamma conditional after the collapsed parameter transitions.

The function does not estimate a mean and does not silently center the observations. Simulation series have mean zero by construction. For applications, the user should choose and document any centering or alternative mean model explicitly.

This target uses zero initial conditions and can be sensitive to initial observations, especially near d=0.5. It should be labelled `BayesConditional`, not `BayesParam` or the original filtered benchmark, in comparisons. No claim is made that it reproduces the manuscript's benchmark values.
