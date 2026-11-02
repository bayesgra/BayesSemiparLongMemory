# Bayesian semiparametric long-memory regression

Research code accompanying Clara Grazian's manuscript **Bayesian Semiparametric
Regression for Dependent Data: An Application to Long Memory**.

This project implements the manuscript's conditional mean-zero, kernel
stick-breaking error mixture for log-periodogram regression. It replaces the
older `longmemory_DP.R` model, whose mixture weights did not depend on frequency
and whose errors were not constrained to have mean zero.

**Status:** a corrected implementation and reproducible workflow, not a
reproduction of the manuscript's numerical results. Full simulations and
substantive chain diagnostics must be rerun before reporting new estimates.
All 42 local tests passed; see the [validation record](docs/validation.md) for
executed workflows and limitations.

## Quick start

Install R, unpack the project, and run these commands from its root directory.
The core code requires no contributed R packages and scripts do not install
anything automatically. The optional Nile application requires `longmemo` for
the original `NileMin` data.

```sh
Rscript tests/run_tests.R
Rscript scripts/run_example.R
Rscript scripts/run_simulations.R --mode=smoke
```

The example and smoke study use deliberately short chains and small truncation
levels to check execution. They are not suitable for scientific inference.
Scripts also work when invoked by absolute path from another directory.
Arguments use `--name=value` syntax.

## Fit your own series

```r
source("R/periodogram.R")
source("R/longmemory_ksb.R")

# Replace series with your finite numeric time series.
dat <- periodogram_data(series, taper = "order1")
fit <- longmemory_ksb_sampler(
  Y = dat$Y, omega = dat$omega,
  H = 30, runs = 30000, burn = 6000, thin = 5,
  seed = 101
)
mean(fit$d)
quantile(fit$d, c(0.025, 0.975))
plot(fit$d, type = "l", ylab = "d", xlab = "Saved draw")
fit$diagnostics
```

Run several independent chains and assess their mixing before using posterior
summaries. The core sampler accepts log-periodogram responses `Y` and increasing
positive Fourier frequencies `omega` in radians. It constructs
`x=log(4*sin(omega/2)^2)` and `t=omega/max(omega)` internally. It estimates both
intercept `c` and memory parameter `d`, with mean `c-d*x`.

The default preprocessing uses the order-one complex taper and retains every
second Fourier frequency within `floor(n^0.8)`. Use
`periodogram_data(series, taper="none")` for the ordinary periodogram. See
[the explicit conventions](docs/reproducibility.md) before comparing analyses.

## Outputs and sampler controls

The fit contains vectors `c`, `d`, `alpha`, `lambda`, and `log_pseudolik`;
matrices `pi`, `mu`, `sigma1_sq`, `sigma2_sq`, `V`, and `Gamma`; the input/design
vectors; settings; and acceptance diagnostics. Rows of parameter matrices are
saved iterations. `terminal_weight_max` tracks the largest remaining-stick weight
among observed covariates at each saved draw. This helps assess sensitivity to
`H`; it is not an automatic convergence test.

`seed=NULL` uses the current R random-number stream. Passing a seed resets that
stream explicitly. Proposal scales can be changed with
`proposal=list(pi=0.35,V=0.4,lambda=0.2)`; these are random-walk scales on logit,
logit and log scales, respectively. No automatic adaptation is performed.

Outputs from scripts are placed under `results/` by default, excluded from git.
Use `--out=/path/to/new_directory` to select a different destination. Use a new
directory for each run; scripts do not resume unfinished jobs automatically.

## Full simulation design

```sh
Rscript scripts/run_simulations.R --mode=paper --out=results/new_paper_run
```

This schedules the manuscript's Gaussian design plus the standardized Student-t
extension: **21,600 fits** at the default 200 replicates per setting. It is a
substantial research computation, not a quick verification command. Defaults in
paper mode are `H=30`, `runs=30000`, `burn=6000`, `thin=5`. Override counts for a
pilot, for example:

```sh
Rscript scripts/run_simulations.R --mode=paper --reps=2 --runs=1000 --burn=500 --thin=1 --out=results/pilot
```

The simulation workflow compares the new model with a clearly labeled base-R
raw-periodogram GPH estimator. That comparator is not certified identical to
`LongMemoryTS::gph`, and the workflow does not reproduce the paper's parametric
Bayesian comparator. Its ARFIMA generator uses finite convolution and warmup;
check warmup sensitivity before drawing conclusions, especially near `d=0.5`.

## Real-data applications

```sh
Rscript scripts/run_real_data.R --dataset=sp500
Rscript scripts/run_real_data.R --dataset=sp500_covid
Rscript scripts/run_real_data.R --dataset=temperature
# Requires the separately installed longmemo package:
Rscript scripts/run_real_data.R --dataset=nile
```

Defaults are the paper priors and `H=30`, with 30,000 iterations. The S&P
workflows calculate log returns from the bundled prices; the subsample selects
returns dated July 30, 2019 through October 31, 2020. The temperature workflow
uses first differences of monthly HadCRUT anomalies, following the old code.
The Nile workflow uses `longmemo::NileMin`, not the different `datasets::Nile`
series. Outputs include the analyzed data, retained frequencies, fit, summary and
diagnostic plots. Application GPH uses the raw grid from index 1.

The Bayesian applications default to tapered/thinned frequencies. Specify
`--taper=none` to use the raw grid instead; the paper does not resolve this
application-versus-simulation convention completely. Neither choice is claimed
to reproduce its old estimates.

After simulations, create performance summaries and plots with:

```sh
Rscript scripts/summarize_simulations.R --input=results/simulation_smoke/estimates.csv
```

An optional [conditional ARFIMA comparator](docs/conditional-benchmark.md) is
also included as `R/arfima_benchmark.R`. It has a coherent conditional likelihood
but is a different benchmark from the paper's cited filtered MCMC method. It is
not included automatically in the main simulation runner.

## Project guide

| Path | Purpose |
|---|---|
| `R/longmemory_ksb.R` | Shared mean-zero kernel stick-breaking sampler and density helpers |
| `R/periodogram.R` | Tapered/raw preprocessing and base-R GPH comparison |
| `R/simulate_arfima.R` | ARFIMA simulation with Gaussian or standardized Student-t innovations |
| `scripts/` | Portable examples, simulations and applications |
| `tests/run_tests.R` | Mathematical identities, input checks and sampler smoke checks |
| `.github/workflows/tests.yml` | Automated R checks on GitHub |
| `data/` | CSV inputs copied from the supplied archive |
| [Model](docs/model.md) | Equations, priors and manuscript mapping |
| [Migration](docs/migration.md) | Changes from the old model and identified workflow bugs |
| [Reproducibility](docs/reproducibility.md) | Conventions, design, diagnostics and replication gaps |

The top-level `longmemory_DP.R` is only a migration loader. Calling the old
`longmemory_dp_sampler()` raises an explanatory error instead of silently
reinterpreting incompatible old arguments. Use `longmemory_ksb_sampler()`.

## Citation and publication

Citation metadata is in `CITATION.cff`. Update it when the manuscript's publication
information becomes available. No code license is assigned here; the author
should select one and check redistribution/attribution conditions for the
bundled data before making the repository public. The supplied manuscript and
old working-session files are not included in this project.
