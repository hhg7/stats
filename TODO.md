# TODO

## Left open

- **A double build keeps the Numerical Recipes fraction below a shape of
  1e3** for the cases R's bratio() sends to bup(), apser() and fpser():
  at a shape of 2^-60 it is 1e16 to 1e18 epsilons out with some logs NaN,
  and up to 2.5e3 in the bup() region. A wide build has taken
  ib_bratio_small() there since 0.3429; doing the same on a double
  changes its results, so it needs R comparisons rather than mpmath's.
  t/incomplete.gamma.beta.mpmath.t skips its "eps" rows on a double.
- **Several grouping factors are one dense matrix** in `lmer` and `glmer`,
  factored whole: once per criterion evaluation in `lmer`, and at every
  PIRLS step in `glmer`, so `grouseticks`' three crossed factors (584
  random effects) took 13 minutes on a desktop (2026-11-22, double NV)
  where `lme4`'s sparse Cholesky takes about a second. The first term (the
  one with the most levels, which `lme4`'s order puts first) has a
  diagonal block of `Z'WZ` when it is scalar, so a Schur complement on it
  would leave a dense factorization of only the other terms' random
  effects (181 for `grouseticks`); a sparse Cholesky with a fill-reducing
  order is the general answer.
  `t/glmer.poisson.R.t`'s `grouse_full` runs only under `EXTENDED_TESTING` for
  this reason.
