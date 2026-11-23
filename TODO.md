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

## Planned

- **`glmer`: generalized linear mixed models**, `lme4::glmer()` and
  `glmer.nb()`, for count and binary outcomes with random effects (wanted
  for a random-intercept negative binomial over child-months, where a
  fixed-effects `glm` drops every child with no events). It would be a
  separate function sharing `lmer`'s formula parser, `LmmTerm`/theta
  mapping and per-level block factorization, with `glm`'s arguments
  (`family`, `offset`, `weights`, `conf_level`, `theta`) and `glm`'s `exp`
  rate ratios. It would not be a `family` argument on `lmer`, which keeps
  its σ, REML and Satterthwaite output unchanged.
  - Method (Bates et al. 2015, J Stat Softw 67(1)): penalized IRLS for the
    conditional modes u at fixed (θ, β), solving
    (Λ′Z′WZΛ + I) u = Λ′Z′W(z − Xβ − offset); Laplace deviance = Σ unit
    deviances (`glm_unit_dev`) + u′u + log|L|². The outer optimization
    follows `lme4`'s two stages: θ alone with β inside PIRLS
    (`nAGQ = 0`), then θ and β jointly. Standard errors come from the
    finite-difference Hessian of the deviance in (θ, β), using the
    Richardson differences `lmer` already has, and the summary reports
    Wald z.
  - What changes in `lmer`'s code: `lmm_prep` forms Z′Z, Z′X and Z′y once,
    so an evaluation never touches the n rows. Under PIRLS the weights
    change every iteration, so Z′WZ is rebuilt from the rows each time
    (still O(n) with one grouping factor). Update the block comment above
    `LmmModel` in the same edit.
  - NB: alternate `nb_theta_ml` with refits at fixed theta, as
    `glmer.nb` does and as `glm` already does for `family => 'negbin'`.
  - Steps, each with an R comparison test like `t/lmer.R.t`:
    1. `offset` and `weights` in `lmer` (it croaks on `offset()` now),
       against `lme4::lmer`.
    2. `glmer`, binomial and Poisson, Laplace, any terms on one grouping
       factor, against `lme4::glmer`.
    3. `family => 'negbin'`, against `lme4::glmer.nb`.
    4. `nAGQ > 1` (adaptive Gauss–Hermite) for a single scalar random
       intercept, which is all `lme4` allows it for, against
       `glmer(nAGQ = 10)` and `GLMMadaptive`.
    5. Several grouping factors, on the dense path with `lmer`'s size
       limit.
