# TODO

## The 2026-10-12 audit

The audit of about 40 XSUBs in `LikeR.xs` on 2026-10-12 found the bugs that
stood in this file. All of them were fixed in 0.3217, the same release that
fixed `roc()`, `auc()` and `auroc()`; `Changes` describes each fix, and
`t/audit.2026-10-12.regressions.t`, `t/bedroc.mpmath.t`, `t/croak.leaks.t`,
`t/hist.R.t`, `t/power_t_test.R.scipy.t` and
`t/friedman_mcnemar_prop_cmh.R.t` pin them.

Three divergences from R were kept on purpose rather than fixed, and README
documents each one:

- **`mcnemar_test()`, paired form**: R requires `nlevels(x) == nlevels(y)`
  after `complete.cases()`; this cross-tabulates over the union of the levels.
- **`sample()`**: `n` defaults to 1, not to a permutation of the whole input,
  and there is no `replace` or `prob`.
- **`ljoin()`, `add_data()`, `transpose()`, `sample()`** copy values, as
  assignment does: a value that is a reference is copied as a reference, so
  the thing it points at is still shared.

## Left open

- **`cmh_test()`** skips a stratum whose total is 0, and runs on one whose
  total is 1. R's `mantelhaen.test()` stops on both ("sample size in each
  stratum must be > 1").
- **`sum()`, `mean()`, `sd()`, `var()`** pay about 0.5 ns per element for the
  compensated sum that matches R's long-double accumulation (0.93 to 1.40
  ns per element for `sum()` on 1e6 NVs). If that matters more than the last
  bits of a badly conditioned sum, `sum()` alone could go back to a plain
  four-lane sum.
- **`rank()`** compares `Stack_off_t` with `size_t` (`opt_start`) in two
  loops, which `-Wextra` reports under strict C99. The default build flags do
  not.
