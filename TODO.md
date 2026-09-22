# TODO

Written 2026-09-22 while checking whether LikeR could run every analysis in
`~/ui/nursing/nursing.proposal.tex` (Table 3 there lists the model for each
aim). The six aims all regress an outcome on hours of home nursing per week.
Each outcome is measured on children who are followed for different lengths of
time and who contribute several years each. What follows is everything that
cannot be done with 0.318 as it stands, roughly in order of how much each item
unblocks.

The R function after each item is the reference to validate against, in the
same way `glm` is already checked against `glm` and `MASS::glm.nb`.

## glm

- **`offset`.** Right now `glm` accepts only `formula`, `data`, `family`,
  `theta` and `conf.level`, and dies on anything else. Count models of
  admissions or ED visits need `offset = log(person-time)` so that a child seen
  for four months is not compared with a child seen for twelve as if they were
  alike. Putting `log(time)` in the formula doesn't do the same job, because
  that estimates its coefficient instead of fixing it at 1. This needs to work
  for `poisson` and `negbin`, and the `negbin` theta search has to see the
  offset too. Reference: `glm(y ~ x + offset(log(t)), family = poisson)` and
  `MASS::glm.nb(y ~ x + offset(log(t)))`.
- **Prior weights (`weights`).** Needed for survey weights in Aim 5 and useful
  generally. Reference: `glm(..., weights = w)`. R warns about non-integer
  weights with `binomial`, and that warning should come through here too.
- **Cluster-robust (sandwich) covariance.** Each child contributes several
  child-years, so the model-based standard errors are too small. The minimum is
  `vcov_type => 'HC0'` and `cluster => 'child_id'`, with the CI, `z` and
  p-value recomputed from the robust covariance. The same machinery gives the
  "modified Poisson" risk ratio for a binary outcome: a `poisson` fit on a 0/1
  response with HC0 standard errors (Zou 2004, Am J Epidemiol 159:702). Aim 5
  uses it. Reference: `sandwich::vcovCL(fit, cluster = ~id, type = "HC0")` and
  `sandwich::vcovHC`, with `lmtest::coeftest` for the tests.
- **Absorbing a high-cardinality factor (fixed effects).** A within-child
  comparison means a factor with thousands of levels. As dummy columns in a
  dense X that is infeasible. For Poisson, conditioning on the child's total is
  identical to including child dummies, so the factor can be absorbed by
  demeaning or iterating within groups rather than expanded. Reference:
  `fixest::fepois(y ~ x | id)`, which should agree with
  `glm(y ~ x + factor(id), family = poisson)` on a small case.
- **Zero-truncated negative binomial, and a hurdle model built on it.** Aim 2
  (inpatient days) is two parts. The first is whether there was any stay,
  which `binomial` already covers. The second is how many days, given at
  least one, which needs a count family truncated at zero. Reference:
  `pscl::hurdle(y ~ x, dist = "negbin", zero.dist = "binomial")`, and
  `countreg::zerotrunc` for the truncated part alone.

## coxph

- **Counting-process (start, stop] input.** A child's hours change during
  follow-up, so hours have to be a time-varying covariate. That needs each
  child split into intervals with their own covariate values. The same input
  handles late entry, since a child enters the cohort on the date nursing is
  authorized and not at birth. Right now `coxph` takes only
  `(\@time, \@status, covariates)`. Reference:
  `survival::coxph(Surv(tstart, tstop, event) ~ x)`.
- **`strata`.** Separate baseline hazards by technology dependence and
  chronic-condition count, as every model in the proposal is stratified.
  Reference: `coxph(... ~ x + strata(g))`.
- **Robust variance with `cluster`.** This is needed as soon as one child
  contributes several intervals. Reference: `coxph(..., cluster = id)`, which
  reports `robust se`.
- **A formula/data interface like `glm`'s.** This is a convenience, not a
  blocker. Once start/stop and strata exist, parallel array refs become hard to
  keep straight.
- **Changepoint profile.** Aim 4 puts a step near 56 hours a week and wants the
  location estimated with a CI. Profiling works: refit with
  `hours > c` for each `c` on a grid and keep `loglik`. That can be done in a
  loop once time-varying covariates exist, so it may not need its own
  function. Reference: `segmented::segmented` or a manual profile in R.

## Models not present at all

- **Linear mixed model** with a random intercept and a random slope per
  subject, fitted by REML, reporting fixed effects with standard errors and
  variance components. Aim 6 (weight-for-age z over time) needs this.
  Reference: `lme4::lmer(z ~ hours * time + (time | id))`, with
  `lmerTest` for the denominator degrees of freedom. The Aim 6 cohort is small,
  so this is the item where falling back to R costs least.
- **Two-stage least squares / instrumental variables.** The instruments in the
  proposal are state rate changes, local nurse supply and waiver waiting lists.
  This needs the second-stage estimate with correct standard errors (not the
  naive standard errors from running two `lm` fits one after the other), the
  first-stage partial F on the excluded instruments, and ideally a control-function
  version for count outcomes. Reference: `AER::ivreg` or `ivreg::ivreg`, with
  `summary(fit, diagnostics = TRUE)` for weak-instrument and Wu-Hausman tests.
- **Nested-model comparison.** `anova` is sequential Type I over the terms of a
  single fit. A joint F (or likelihood-ratio) test of a block of terms needs
  `anova(fit_small, fit_big)`, which is also how the first-stage F would be
  computed before IV exists. Reference: `anova(m0, m1)` for `lm`, and
  `anova(m0, m1, test = "Chisq")` for `glm`.

## Survey design

- **Design-based estimation for weighted surveys:** weights, strata and
  clusters (PSUs), with standard errors that respect them. Aim 5's survey
  will carry nonresponse weights, and the National Survey of Children's Health,
  used as its benchmark, is published with design variables. `weights` on
  `glm` alone gives correct point estimates but wrong standard errors.
  Reference: `survey::svydesign` with `survey::svyglm`.

## Not LikeR's to fix, but it decides how much of the above matters

The claims aims run inside the CMS Virtual Research Data Center, where
installing outside compiled software may not be allowed. Ask ResDAC before
building anything above only for the T-MSIS work. If XS modules can't go in,
these items still matter for the HCUP analyses and anything else done outside
the enclave.

## Status, 2026-09-22 (0.319)

Every item above except the last section is in 0.319; see `Changes` for what
each became and which `t/*.R.t` file pins it to the R reference. Two were done
as documentation rather than new functions, as the items themselves allowed:
the changepoint profile is a loop over `coxph` fits (README.md, under
`coxph`), and the control-function IV for count outcomes is a first-stage `lm`
plus a `glm` on its residual (README.md, under `glm`), whose second-stage
standard errors need a bootstrap. `svyglm` handles one-stage designs only.
