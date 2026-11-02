# Validation record — 2 November 2026

Validation used R 4.3.3 on Linux. The final test suite passed all 42 checks,
including numerical integration of the mean-zero atom, weight normalization,
the transformed Beta prior target, FFT normalization, the regression sign,
seed reproducibility, sampler support and conditional ARFIMA filtering.

Short end-to-end runs completed for the example, Gaussian and Student-t
simulations, performance summaries, full S&P 500, its COVID-period subsample,
and monthly temperature differences. Their CSV, RDS and PDF outputs were
generated successfully. The Nile runner's missing-package error was verified;
its full workflow was not executed because `longmemo` was unavailable.
GitHub Actions is configured for three operating systems, but those remote
jobs have not yet run.

Four independent 2,500-iteration chains were additionally run on one synthetic
asymmetric mean-zero regression dataset with 160 observations and true d=0.25.
With H=3, 1,000 burn-in iterations and thinning 2, all saved values were finite.
The four posterior means for d were 0.3355, 0.3313, 0.3326 and 0.3314; all four
95% intervals contained 0.25. This single dataset cannot assess bias or coverage.
The deliberately small truncation retained substantial terminal mass, and the
nuisance parameters require further mixing and truncation checks.

Independent code review checked the allocation, regression, linked atom mean,
variance, stick, knot, concentration and bandwidth updates. A Beta-prior
Jacobian error found during development was corrected and covered by a
regression test before this validation.

These checks establish implementation and execution checks, not convergence,
frequentist calibration under dependent observations, or reproduction of the
paper's tables. The full 21,600-fit study has not been run.

## Application count discrepancy

For the supplied S&P data, the documented cutoff `floor(n^0.8)` produces 301
raw ordinates and 150 tapered/thinned ordinates. The manuscript reports 302;
the old code used `floor(1+n^0.8)`. The new project follows the manuscript's
formula and records actual retained indices with every fit.
