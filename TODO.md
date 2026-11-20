# TODO

## Left open

- **`pnorm()` and `qnorm()` on a long-double or `__float128` perl** are
  still double-accurate. 0.3426 made the incomplete gamma and beta compute
  at the NV's width there (`LIKER_WIDE_NV` in `LikeR.xs`), but `c_pnorm()`
  and `approx_pnorm()` are double algorithms with a cut-off at z = -38.47,
  and `M_PI` is still the double literal outside the paths 0.3426 changed.
  The wide gamma path has its own `pg_pnorm()`/`pg_dpnorm()` (erfcx-based),
  which a wide `pnorm()` could build on.
- **The wide incomplete gamma and beta are short of the epsilon in three
  places** (0.3426's `Changes` entry has the measurements): 25 to 39
  epsilons of a `__float128` for a < 1 with x between 1 and 3 (the
  continued fraction's 400 terms of rounding), 86 for a shape of 3/2
  against 1e4 and up (`bgrat()`'s expansion), and TOMS 708's `bup()`,
  `apser()` and `fpser()` are not ported, so a maximum shape of at most 15
  with x0 b0 > 1.5 and x0 in [0.1, 0.29) stays on the Numerical Recipes
  fraction. Nothing past a shape of 1e7 (beta) or 1e10 (gamma) has an
  mpmath reference; the Temme path at 1e11 to 1e15 is covered only by the
  R-derived tests.
