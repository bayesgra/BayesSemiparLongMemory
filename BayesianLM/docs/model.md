# Model and manuscript mapping

This implementation follows Sections 2, 4 and 5 of the supplied manuscript,
*Bayesian Semiparametric Regression for Dependent Data: An Application to Long
Memory*, by Clara Grazian. It replaces the old, covariate-independent Dirichlet
process mixture. It does not certify or reproduce the manuscript's numerical
results.

## Regression and inferential target

For retained positive Fourier frequencies `omega[j]`,

```text
x[j] = log(4 * sin(omega[j]/2)^2)
t[j] = omega[j] / max(omega)
Y[j] = c - d*x[j] + U[j]
u[j] = Y[j] - c + d*x[j]
```

`Y` contains the log-periodogram, not the periodogram itself. The intercept `c`
absorbs the constant centering of log-periodogram errors. The unknown
short-memory bias in manuscript equation (10) is omitted from the fitted model,
as specified in Section 4. The frequency covariate is normalized by its largest
retained value, not min-max scaled.

The target is proportional to

```text
prior(parameters) * product_j f(Y[j] - c + d*x[j] | t[j]).
```

This is the marginal **pseudo-posterior**. The code does not specify a joint
probability model for dependence among periodogram ordinates. There is no extra
`1/m` power on the product. Posterior quantiles are pseudo-posterior credible
intervals; they are not automatically calibrated confidence intervals under
dependence.

## Mean-zero errors: equations (2) and (3)

For each atom `h`,

```text
b_h(u) = pi[h] * Normal(u; mu[h], variance1[h])
       + (1-pi[h]) * Normal(u; -pi[h]*mu[h]/(1-pi[h]), variance2[h])
f(u | t) = sum_h W_h(t) * b_h(u).
```

Each atom has mean exactly zero because its two weighted means cancel. Thus
frequency-dependent weights cannot change the regression mean. This is stronger
than assigning a mean-zero prior to unconstrained component locations. Each atom
has two distinct variances; with `H=30`, the density can contain 60 normal terms.

For `h < H`, define

```text
Z_h(t) = V[h] * K(t, Gamma[h], lambda)
W_h(t) = Z_h(t) * product_{l<h}(1-Z_l(t)).
```

The final weight is the entire remaining stick,

```text
W_H(t) = product_{l<H}(1-Z_l(t)).
```

This is the finite approximation specified in Section 5. The final atom has no
ordinary stick or knot. Simply setting `V[H]=1` while multiplying it by a kernel
would leave missing probability mass. This code does not sample the infinite
mixture of Section 2 exactly.

## Default priors from Section 5

| Quantity | Prior or default |
|---|---|
| Number of atoms `H` | 30 |
| `c` | Normal(0, variance 100), restricted to [-50,50] |
| `d` | Uniform(-1, 0.5) |
| `V[h]`, `h < H` | Beta(1, alpha) |
| `alpha` | Gamma(shape 2, rate 1) |
| `Gamma[h]` | Uniform(0,1) |
| Shared `lambda` | Gamma(shape 2, rate 1) |
| `pi[h]` | Beta(2,2) |
| `mu[h]` | Normal(0, variance 10), restricted to [-10,10] |
| Both atom variances | Inverse-Gamma(shape 2, scale 2), restricted to [0.1,50] |

Here inverse-Gamma density is proportional to
`v^(-a-1)*exp(-b/v)`, equivalently `1/rgamma(..., shape=a, rate=b)`.
The variance bounds apply to variance, not standard deviation. The priors are
independent apart from the specified dependence of the sticks on `alpha`.

## Interpretation limits

The theoretical results impose additional support, separation and concentration
conditions. Running this code does not establish those conditions. In particular,
the computational Gamma prior on the bandwidth is not bounded away from zero,
and `pi` arbitrarily close to one can produce a large second normal mean even
when `mu` is bounded. The manuscript's theoretical bounds should not be inferred
from the computational defaults alone.

Inspect the final component's weight across the covariate domain and rerun with
larger `H`; appreciable terminal mass signals that truncation sensitivity deserves
attention. A small terminal mass is not, by itself, a proof of convergence or
adequate posterior exploration.

## Posterior computation

The sampler augments each observation with its outer mixture allocation and its
choice of the two normals within that atom. Their joint categorical draw uses
log probabilities. Given these allocations:

- `c` and `d` have a joint Gaussian conditional, restricted to their prior
  rectangle. Up to 1,000 independent Gaussian proposals are tried; if none lies
  inside, the old pair is retained. The resulting state-independent holding
  probability preserves the conditional target, but poor acceptance requires
  investigation.
- `mu` has a truncated Gaussian conditional. Each variance has a truncated
  inverse-Gamma conditional.
- `pi` and `V` use Gaussian random walks on their logit scales, with the
  transformation Jacobians included in Metropolis ratios.
- Knots use independent Uniform(0,1) Metropolis proposals. Shared `lambda` uses
  a Gaussian random walk on the log scale with its Jacobian included.
- `alpha` is sampled from
  `Gamma(shape=2+H-1, rate=1-sum(log(1-V)))`.

All steps target the same finite-mixture pseudo-posterior. The dependence
interpretation does not change these computational updates. Proposal tuning and
run lengths are implementation choices, not additional methodology specified by
the supplied paper.
