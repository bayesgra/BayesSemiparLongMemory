# Migration from the supplied code archive

The supplied `longmemory_DP.R` and its copies in the real-data directories do not
implement the current manuscript. Reusing their fitted draws with the new model
is not valid.

| Old implementation | Current requirement |
|---|---|
| One set of weights for all frequencies | Kernel stick-breaking weights varying with normalized frequency |
| Unconstrained Gaussian atom locations | Two-Gaussian atoms with mean exactly zero |
| One shared residual variance | Two variances per atom |
| No explicit regression intercept | Infer `c` separately from `d` |
| Mean `-2*d*log(abs(omega))` | Mean `c-d*log(4*sin(omega/2)^2)` |
| Unclosed finite stick | Assign all remaining mass to atom `H` |
| Fixed stick hyperparameters | Infer `alpha` under Gamma(2,1) |
| No knots or bandwidth | Infer knots and shared kernel bandwidth |
| Unconditional internal seed reset | Explicit optional seed, controlled by the caller |

The replacement is not a drop-in numerical continuation of old chains. In
particular, old fields such as `xi` and a single `sigma2` no longer describe the
same statistical objects. Import the new model directly and start new chains.

## Other problems identified in the old workflows

- The simulation driver contains a `sourve` typo, missing or incorrect source
  calls, an undefined `d_vec`, and a call passing already logged responses to a
  function that logs its input again.
- The temperature workflow references undefined `ynew` and `omeganew` objects.
- The full S&P workflow can calculate `fft(NileMin)` while using the financial
  sample length. The subsample workflow references missing code and can reuse the
  full-sample cutoff.
- The old filtered ARFIMA benchmark retains target-density caches after changing
  conditioning parameters. This requires a separate correction and validation.

New scripts build their inputs explicitly and use the shared implementation.
The original archive remains the provenance source; working-session files and
old output objects are not inputs to the new analysis.

## Existing results

None of the manuscript's tables, point estimates or coverage claims should be
replaced merely because the source code has been corrected. Run the corrected
model, assess chains and truncation sensitivity, and regenerate the results.
Differences from the earlier output are expected because the statistical model
has changed. See `reproducibility.md` for the distinction between demonstration
runs and full replication.
