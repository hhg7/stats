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

## The 2026-10-20 audit

Six auditors read the XSUBs that no earlier audit had covered. The survival
functions were fixed first; everything else they found was fixed in the same
release, 0.331, and `Changes` describes each fix. The `t/audit.2026-10-20.*.t`
files, `t/distributions.limits.R.t`, `t/scale.R.t`, `t/prcomp.R.t` and the
extended `t/<fn>.R.t` files pin them.

Divergences from R that were kept on purpose, each documented in README or in
its test file:

- **`ks_test()`**: the two-sample one-sided asymptotic p is `exp(-2nq^2)`
  computed directly, where R's `1 - (1 - exp())` rounds anything below about
  1e-16 to 0. Where R's Birnbaum-Tingey sum takes the log of a negative
  number and returns NaN, `ks_test` drops that term and returns 1.
- **`fisher_test()`**: an r x c p-value can differ from R's by up to about
  3e-12, because FEXACT merges paths whose weights agree to 3.45e-7; this one
  is closer to the exact rational value. An all-zero r x c table still croaks,
  as R does, where SciPy gives 1.
- **`zerotrunc()` / `hurdle()`**: a negbin fit whose likelihood rises to the
  Poisson limit reports `theta` as Inf, where countreg's BFGS stops near 1e5.
- **`cmh_test()`**: K = 1 is accepted, where R refuses it.
- **`p_adjust()`**: undef is a p-value of 1 counted in the family, as README
  documents; R would treat it as NA.

## Left open

- **`binom_test(3, 1e18)`** runs an O(n) loop and effectively never returns.
  R fails there too, trying to allocate 37 GB.
- **The incomplete gamma and beta** (0.331's ports of R's pgamma.c and
  TOMS 708) use R's double-precision series constants, so `pchisq`, `pt`,
  `pf` and every p-value built on them are about double-accurate on
  long-double and quadmath perls, not accurate to the NV's width.
- **`fisher_test()`** refuses PR#4688 (4x3, N = 16442), as R does, after
  0.6 s. A table whose total passes 2^31 - 1 croaks where a C long is 32 bits.
- **`sum()`, `mean()`, `sd()`, `var()`** pay about 0.5 ns per element for the
  compensated sum that matches R's long-double accumulation (0.93 to 1.40
  ns per element for `sum()` on 1e6 NVs). If that matters more than the last
  bits of a badly conditioned sum, `sum()` alone could go back to a plain
  four-lane sum.
- **Formulas** (every function that reads one through `lm_formula_terms()`)
  have no parentheses: `y ~ (a + b)*c` and `y ~ (a + b + c)^2` die with "0
  degrees of freedom", as does arithmetic inside a function or `I()`, such as
  `log(x + 1)` or `I(a - b)`. `car`'s own `vif()` example is
  `(income + education)*type`, which `t/vif.R.t` has to write out as its
  expansion. Term removal (`- b`) works since 0.3218.
- **`vif()`** has only `car::vif()`'s default `type = "terms"`; `type =
  "predictor"`, which `car` suggests for a model with interactions, is not
  implemented.
- **`read_table()` of an `.xlsx`** returns its text as undecoded UTF-8 bytes:
  sheet names, the keys of the hash it returns for a workbook, and string
  cells alike. A workbook that `write_table()` wrote from the keys `café` and
  `αβ` reads back with keys of 5 and 4 bytes (`caf\xC3\xA9`, `\xCE\xB1\xCE\xB2`),
  so `$book->{"\x{3b1}\x{3b2}"}` finds nothing, and a cell `α` reads back as 2
  bytes. README's "the two round-trip" for a workbook holds only for ASCII
  keys. `write_table()` itself writes the names correctly (openpyxl 3.1.5
  reads them back as `café` and `αβ`); the fix belongs in `read_table()`,
  which would have to decode every name and string cell, and that changes
  what it returns for every non-ASCII `.xlsx` cell, not only for sheet names.
  Found on 2026-10-26, checking 0.341's workbooks.
