#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Low-level math helpers                                              */
/* ------------------------------------------------------------------ */

/* log( C(n,k) ) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)          */
static double
lchoose(double n, double k)
{
    if (k < 0 || k > n) return -INFINITY;
    if (k == 0 || k == n) return 0.0;
    return lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
}

/* log P(X = x) under Hypergeometric(m, n, k)
 * x balls drawn from m white + n black, k draws total.
 * Matches R: dhyper(x, m, n, k, log=TRUE)                              */
static double
dhyper_log(int x, int m, int n, int k)
{
    if (x < 0 || x > m || x > k || k - x > n) return -INFINITY;
    return lchoose(m, x) + lchoose(n, k - x) - lchoose(m + n, k);
}

/* ------------------------------------------------------------------ */
/*  Noncentral hypergeometric helpers (used for OR estimation / CI)    */
/*  Matches R's inner functions in fisher.test:                         */
/*    dnhyper, mnhyper, pnhyper                                         */
/* ------------------------------------------------------------------ */

/* Compute unnormalised log-weights w_i = logdc[i] + i*log(ncp),
 * then normalise to proper probabilities.
 * 'support' is lo..hi (integers), len = hi - lo + 1.
 * Result stored into 'out' (caller-allocated, length len).            */
static void
dnhyper(const int *support, const double *logdc, int len,
        double ncp, double *out)
{
    int i;
    /* ncp=0: all mass sits on x=lo; avoids 0*(-Inf)=NaN */
    if (ncp == 0.0) {
        for (i = 0; i < len; i++) out[i] = (i == 0) ? 1.0 : 0.0;
        return;
    }
    double log_ncp = log(ncp);
    double maxval  = -INFINITY;

    /* log-sum-exp trick for numerical stability */
    for (i = 0; i < len; i++) {
        out[i] = logdc[i] + support[i] * log_ncp;
        if (out[i] > maxval) maxval = out[i];
    }

    double sum = 0.0;
    for (i = 0; i < len; i++) {
        out[i] = exp(out[i] - maxval);
        sum += out[i];
    }
    for (i = 0; i < len; i++)
        out[i] /= sum;
}

/* E[X] under noncentral hypergeometric with ncp.
 * Matches R: mnhyper <- function(ncp) sum(support * dnhyper(ncp))     */
static double
mnhyper(const int *support, const double *logdc, int len, double ncp)
{
    double *d = (double *)malloc(len * sizeof(double));
    dnhyper(support, logdc, len, ncp, d);
    double s = 0.0;
    int i;
    for (i = 0; i < len; i++)
        s += support[i] * d[i];
    free(d);
    return s;
}

/* P(X <= q)  or  P(X >= q)  under noncentral hypergeometric.
 * upper_tail=0 => P(X <= q),  upper_tail=1 => P(X >= q).
 * Matches R: pnhyper <- function(q, ncp=1, upper.tail=FALSE)           */
static double
pnhyper(int q, const int *support, const double *logdc, int len,
        double ncp, int upper_tail)
{
    double *d = (double *)malloc(len * sizeof(double));
    dnhyper(support, logdc, len, ncp, d);
    double s = 0.0;
    int i;
    for (i = 0; i < len; i++) {
        if (upper_tail) {
            if (support[i] >= q) s += d[i];
        } else {
            if (support[i] <= q) s += d[i];
        }
    }
    free(d);
    return s;
}

/* ------------------------------------------------------------------ */
/*  Brent's method root-finder  (R's uniroot equivalent)               */
/* ------------------------------------------------------------------ */

typedef struct {
    const int    *support;
    const double *logdc;
    int           len;
    double        target;   /* RHS of equation f(ncp) - target = 0     */
    int           use_mnhyper;   /* 1 => mnhyper,  0 => pnhyper variant */
    int           q;             /* for pnhyper variants                 */
    int           upper_tail;
} BrentCtx;

static double
brent_f(double ncp, BrentCtx *ctx)
{
    if (ctx->use_mnhyper)
        return mnhyper(ctx->support, ctx->logdc, ctx->len, ncp) - ctx->target;
    else
        return pnhyper(ctx->q, ctx->support, ctx->logdc, ctx->len,
                       ncp, ctx->upper_tail) - ctx->target;
}

/* Returns NAN if no root found within [a,b] with same sign convention. */
static double
brent_root(double a, double b, double tol, BrentCtx *ctx)
{
    double fa = brent_f(a, ctx);
    double fb = brent_f(b, ctx);

    /* If both endpoints have the same sign, clamp to the nearer bound. */
    if (fa * fb > 0.0) {
        return (fabs(fa) < fabs(fb)) ? a : b;
    }

    double c = a, fc = fa, d = b - a, e = d;
    int iter;

    for (iter = 0; iter < 1000; iter++) {
        if (fb * fc > 0.0) { c = a; fc = fa; d = e = b - a; }
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        double tol1 = 2.0 * DBL_EPSILON * fabs(b) + 0.5 * tol;
        double xm   = 0.5 * (c - b);
        if (fabs(xm) <= tol1 || fb == 0.0) return b;

        double p, q_val, r;
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            double s = fb / fa;
            if (a == c) {
                p = 2.0 * xm * s;
                q_val = 1.0 - s;
            } else {
                q_val = fa / fc;
                r     = fb / fc;
                p = s * (2.0 * xm * q_val * (q_val - r) - (b - a) * (r - 1.0));
                q_val = (q_val - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q_val = -q_val; else p = -p;
            if (2.0 * p < 3.0 * xm * q_val - fabs(tol1 * q_val) &&
                2.0 * p < fabs(e * q_val)) {
                e = d; d = p / q_val;
            } else {
                d = xm; e = d;
            }
        } else {
            d = xm; e = d;
        }
        a = b; fa = fb;
        if (fabs(d) > tol1) b += d; else b += (xm > 0) ? tol1 : -tol1;
        fb = brent_f(b, ctx);
    }
    return b;
}

/* ------------------------------------------------------------------ */
/*  MLE of odds ratio (ncp) for 2x2 table                              */
/*  Matches R: mle <- function(x) { ... uniroot(mnhyper - x) ... }    */
/* ------------------------------------------------------------------ */
static double
mle_ncp(const int *support, const double *logdc, int len, int observed)
{
    int lo = support[0];
    int hi = support[len - 1];

    if (observed == lo) return 0.0;
    if (observed == hi) return INFINITY;

    BrentCtx ctx;
    ctx.support      = support;
    ctx.logdc        = logdc;
    ctx.len          = len;
    ctx.target       = (double)observed;
    ctx.use_mnhyper  = 1;
    ctx.q            = 0;
    ctx.upper_tail   = 0;

    /* start at 1e-10 to avoid the ncp=0 boundary */
    return brent_root(1e-10, 1e7, 1e-7, &ctx);
}

/* ------------------------------------------------------------------ */
/*  Confidence interval for odds ratio                                  */
/*  Matches R's ci.ncp inner function in fisher.test                   */
/* ------------------------------------------------------------------ */
static void
ci_ncp(const int *support, const double *logdc, int len,
       int observed, double conf_level, const char *alternative,
       double *ci_lo, double *ci_hi)
{
    int lo = support[0];
    int hi = support[len - 1];
    double alpha = 1.0 - conf_level;

    BrentCtx ctx;
    ctx.support     = support;
    ctx.logdc       = logdc;
    ctx.len         = len;
    ctx.use_mnhyper = 0;

    if (strcmp(alternative, "less") == 0) {
        *ci_lo = 0.0;
        /* upper: pnhyper(observed, ncp, upper.tail=FALSE) == alpha */
        if (observed == hi) {
            *ci_hi = INFINITY;
        } else {
            ctx.q          = observed;
            ctx.upper_tail = 0;
            ctx.target     = alpha;
            *ci_hi = brent_root(1e-10, 1e7, 1e-7, &ctx);
        }
    } else if (strcmp(alternative, "greater") == 0) {
        *ci_hi = INFINITY;
        /* lower: pnhyper(observed, ncp, upper.tail=TRUE) == alpha */
        if (observed == lo) {
            *ci_lo = 0.0;
        } else {
            ctx.q          = observed;
            ctx.upper_tail = 1;
            ctx.target     = alpha;
            *ci_lo = brent_root(1e-10, 1e7, 1e-7, &ctx);
        }
    } else {
        /* two.sided: alpha/2 each tail */
        double a2 = alpha / 2.0;

        if (observed == lo) {
            *ci_lo = 0.0;
        } else {
            ctx.q          = observed;
            ctx.upper_tail = 1;
            ctx.target     = a2;
            *ci_lo = brent_root(1e-10, 1e7, 1e-7, &ctx);
        }

        if (observed == hi) {
            *ci_hi = INFINITY;
        } else {
            ctx.q          = observed - 1;    /* pnhyper(x-1, ncp, upper.tail=FALSE) */
            ctx.upper_tail = 0;
            ctx.target     = a2;
            *ci_hi = brent_root(1e-10, 1e7, 1e-7, &ctx);
        }
    }
}

/* ------------------------------------------------------------------ */
/*  P-value computation                                                 */
/*  Matches R's internal p-value logic for 2x2 fisher.test             */
/* ------------------------------------------------------------------ */

/* For "less":    P(X <= x_obs)   under H0 (ncp=1)                     */
/* For "greater": P(X >= x_obs)                                         */
/* For "two.sided": sum of dhyper(x) for all x where                   */
/*                  dhyper(x) <= dhyper(x_obs)  (R's default)          */
static double
fisher_pvalue(const double *logdc, int len, int obs_idx,
              const char *alternative)
{
    double pval;
    int i;

    if (strcmp(alternative, "less") == 0) {
        pval = 0.0;
        for (i = 0; i <= obs_idx; i++)
            pval += exp(logdc[i]);
    } else if (strcmp(alternative, "greater") == 0) {
        pval = 0.0;
        for (i = obs_idx; i < len; i++)
            pval += exp(logdc[i]);
    } else {
        /* two.sided: sum where log-prob <= observed log-prob + tiny tol */
        double obs_lp = logdc[obs_idx];
        pval = 0.0;
        for (i = 0; i < len; i++) {
            if (logdc[i] <= obs_lp + 1e-7)
                pval += exp(logdc[i]);
        }
    }
    if (pval > 1.0) pval = 1.0;
    return pval;
}

/* ------------------------------------------------------------------ */
/*  Main 2x2 Fisher exact test                                          */
/* ------------------------------------------------------------------ */

/* Mirrors R's fisher.test for a 2x2 matrix.
 * Returns a hash-ref with: p_value, odds_ratio, conf_int (arrayref),
 * alternative, method, statistic (the table cell [0,0]).              */
static HV *
fisher_test_2x2(int a, int b, int c, int d,
                const char *alternative, double conf_level,
                int do_conf_int)
{
    /* Marginals  (R notation: m = row1 total, n = row2 total,
     *             k = col1 total;  support of X = a)                  */
    int m  = a + b;    /* row 1 total  */
    int n  = c + d;    /* row 2 total  */
    int k  = a + c;    /* col 1 total  */
    int lo = (k - n > 0) ? k - n : 0;
    int hi = (k < m)     ? k     : m;
    int len = hi - lo + 1;

    /* Build support array and log-probability table */
    int    *support = (int *)    malloc(len * sizeof(int));
    double *logdc   = (double *) malloc(len * sizeof(double));
    int i;
    for (i = 0; i < len; i++) {
        support[i] = lo + i;
        logdc[i]   = dhyper_log(lo + i, m, n, k);
    }
    /* Index of observed cell value in support */
    int obs_idx = a - lo;
    /* P-value */
    double pval = fisher_pvalue(logdc, len, obs_idx, alternative);
    /* Odds ratio (MLE) */
    double or_est = mle_ncp(support, logdc, len, a);
    /* Confidence interval */
    double ci_lo = 0.0, ci_hi = INFINITY;
    if (do_conf_int) {
        ci_ncp(support, logdc, len, a, conf_level, alternative,
               &ci_lo, &ci_hi);
    }
    free(support);
    free(logdc);
    /* Build return hash */
    HV *result = newHV();
    hv_store(result, "p_value",      7, newSVnv(pval),    0);
    hv_store(result, "odds_ratio",   10, newSVnv(or_est), 0);
    hv_store(result, "alternative",  11, newSVpv(alternative, 0), 0);
    hv_store(result, "method",       6,
             newSVpv("Fisher's Exact Test for Count Data", 0), 0);
    hv_store(result, "statistic",    9, newSViv(a), 0);
    hv_store(result, "conf_level",   10, newSVnv(conf_level), 0);

    /* conf.int as [ lo, hi ] */
    if (do_conf_int) {
        AV *ci = newAV();
        av_push(ci, newSVnv(ci_lo));
        av_push(ci, (ci_hi == INFINITY) ? newSVpv("Inf", 0)
                                        : newSVnv(ci_hi));
        hv_store(result, "conf_int", 8, newRV_noinc((SV *)ci), 0);
    }

    return result;
}

/* ------------------------------------------------------------------ */
/*  XS interface                                                        */
/* ------------------------------------------------------------------ */

MODULE = Statistics::FisherTest    PACKAGE = Statistics::FisherTest

PROTOTYPES: DISABLE

# fisher_test(a, b, c, d, alternative, conf_level, do_conf_int)
#
# Arguments:
#   a, b, c, d   - the four cells of a 2x2 contingency table:
#                    | a  b |
#                    | c  d |
#   alternative  - "two.sided", "less", or "greater"
#   conf_level   - confidence level for CI, e.g. 0.95
#   do_conf_int  - 1 to compute CI, 0 to skip

SV *
fisher_test(a, b, c, d, alternative, conf_level, do_conf_int)
    int    a
    int    b
    int    c
    int    d
    char  *alternative
    double conf_level
    int    do_conf_int
  CODE:
    HV *result = fisher_test_2x2(a, b, c, d,
                                  alternative, conf_level, do_conf_int);
    RETVAL = newRV_noinc((SV *)result);
  OUTPUT:
    RETVAL
