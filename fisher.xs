/* ===========================================================================
 * Stats::LikeR  --  Fisher's Exact Test (2x2), R-compatible
 *
 * This file fixes two categories of bug:
 *   (1) the statistics: the conditional-MLE odds ratio and the
 *       noncentral-hypergeometric confidence interval, matching R's
 *       stats::fisher.test (including R's uniroot/zeroin tolerance, which
 *       is what makes the published CI bounds reproducible).
 *   (2) the Perl<->C glue: signed integers, input validation, and fetching
 *       each hash value exactly once.
 *
 * Drop the C helper block into your .xs (above MODULE) and replace the
 * fisher_test body with the one below.
 * =========================================================================== */

#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ helpers */

#define FT_EPS 2.220446049250313e-16
#define FT_TOL 0.0001220703125  /* .Machine$double.eps^0.25, R uniroot default */

static double ft_lchoose(long n, long k) {
    if (k < 0 || k > n || n < 0) return -INFINITY;
    return lgamma((double)n + 1) - lgamma((double)k + 1) - lgamma((double)(n - k) + 1);
}

typedef struct {
    long lo, hi, ns, m, n, k, x;
    double *logdc;   /* central log hypergeometric density over the support */
} ft_support;

static int ft_init(ft_support *S, long a, long b, long c, long d) {
    S->m = a + c; S->n = b + d; S->k = a + b; S->x = a;
    S->lo = (S->k - S->n > 0) ? (S->k - S->n) : 0;
    S->hi = (S->k < S->m) ? S->k : S->m;
    S->ns = S->hi - S->lo + 1;
    if (S->ns <= 0) { S->logdc = NULL; return 0; }
    Newx(S->logdc, S->ns, double);
    for (long i = 0; i < S->ns; i++) {
        long j = S->lo + i;
        S->logdc[i] = ft_lchoose(S->m, j) + ft_lchoose(S->n, S->k - j)
                    - ft_lchoose(S->m + S->n, S->k);
    }
    return 1;
}
static void ft_free(ft_support *S) { Safefree(S->logdc); S->logdc = NULL; }

static void ft_dnhyper(const ft_support *S, double ncp, double *out) {
    double lncp = log(ncp), mx = -INFINITY;
    for (long i = 0; i < S->ns; i++) {
        out[i] = S->logdc[i] + lncp * (double)(S->lo + i);
        if (out[i] > mx) mx = out[i];
    }
    double s = 0;
    for (long i = 0; i < S->ns; i++) { out[i] = exp(out[i] - mx); s += out[i]; }
    for (long i = 0; i < S->ns; i++) out[i] /= s;
}

static double ft_mnhyper(const ft_support *S, double ncp, double *scratch) {
    if (ncp == 0)     return (double)S->lo;
    if (isinf(ncp))   return (double)S->hi;
    ft_dnhyper(S, ncp, scratch);
    double mu = 0;
    for (long i = 0; i < S->ns; i++) mu += (double)(S->lo + i) * scratch[i];
    return mu;
}

/* upper != 0 => P(X >= q), upper == 0 => P(X <= q) */
static double ft_pnhyper(const ft_support *S, long q, double ncp, int upper, double *scratch) {
    if (ncp == 1.0) {
        double s = 0;
        for (long i = 0; i < S->ns; i++) {
            long j = S->lo + i;
            if (upper ? (j >= q) : (j <= q)) s += exp(S->logdc[i]);
        }
        return s;
    }
    if (ncp == 0.0)   return upper ? (double)(q <= S->lo) : (double)(q >= S->lo);
    if (isinf(ncp))   return upper ? (double)(q <= S->hi) : (double)(q >= S->hi);
    ft_dnhyper(S, ncp, scratch);
    double s = 0;
    for (long i = 0; i < S->ns; i++) {
        long j = S->lo + i;
        if (upper ? (j >= q) : (j <= q)) s += scratch[i];
    }
    return s;
}

/* R's src/library/stats/src/zeroin.c (Brent-Dekker) */
typedef double (*ft_fn)(double t, void *ctx);
static double ft_zeroin(double ax, double bx, ft_fn f, void *ctx, double tol, int maxit) {
    double a = ax, b = bx, fa = f(a, ctx), fb = f(b, ctx), c = a, fc = fa;
    while (maxit-- > 0) {
        double prev = b - a;
        if (fabs(fc) < fabs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
        double tol_act = 2 * FT_EPS * fabs(b) + tol / 2;
        double step = (c - b) / 2;
        if (fabs(step) <= tol_act || fb == 0.0) return b;
        if (fabs(prev) >= tol_act && fabs(fa) > fabs(fb)) {
            double cb = c - b, p, q;
            if (a == c) { double t1 = fb / fa; p = cb * t1; q = 1.0 - t1; }
            else {
                double q0 = fa / fc, t1 = fb / fc, t2 = fb / fa;
                p = t2 * (cb * q0 * (q0 - t1) - (b - a) * (t1 - 1.0));
                q = (q0 - 1.0) * (t1 - 1.0) * (t2 - 1.0);
            }
            if (p > 0) q = -q; else p = -p;
            if (p < 0.75 * cb * q - fabs(tol_act * q) / 2 && p < fabs(prev * q / 2)) step = p / q;
        }
        if (fabs(step) < tol_act) step = step > 0 ? tol_act : -tol_act;
        a = b; fa = fb; b += step; fb = f(b, ctx);
        if ((fb > 0) == (fc > 0)) { c = a; fc = fa; }
    }
    return b;
}

typedef struct { const ft_support *S; double target; double *scratch; int mode; } ft_rc;
/* mode 0: mnhyper(t)-target      1: mnhyper(1/t)-target
   mode 2: pnhyper(x,t,low)-tgt   3: pnhyper(x,1/t,low)-tgt
   mode 4: pnhyper(x,t,up)-tgt    5: pnhyper(x,1/t,up)-tgt */
static double ft_rootf(double t, void *ctx) {
    ft_rc *r = (ft_rc *)ctx; const ft_support *S = r->S;
    switch (r->mode) {
        case 0: return ft_mnhyper(S, t, r->scratch) - r->target;
        case 1: return ft_mnhyper(S, 1.0 / t, r->scratch) - r->target;
        case 2: return ft_pnhyper(S, S->x, t, 0, r->scratch) - r->target;
        case 3: return ft_pnhyper(S, S->x, 1.0 / t, 0, r->scratch) - r->target;
        case 4: return ft_pnhyper(S, S->x, t, 1, r->scratch) - r->target;
        default:return ft_pnhyper(S, S->x, 1.0 / t, 1, r->scratch) - r->target;
    }
}

static double exact_p_value(long a, long b, long c, long d, const char *alt) {
    ft_support S;
    if (!ft_init(&S, a, b, c, d)) return 1.0;
    double *sc; Newx(sc, S.ns, double);
    double p;
    if (!strcmp(alt, "less"))         p = ft_pnhyper(&S, S.x, 1.0, 0, sc);
    else if (!strcmp(alt, "greater")) p = ft_pnhyper(&S, S.x, 1.0, 1, sc);
    else {
        ft_dnhyper(&S, 1.0, sc);
        double dx = sc[S.x - S.lo], relErr = 1 + 1e-7, s = 0;
        for (long i = 0; i < S.ns; i++) if (sc[i] <= dx * relErr) s += sc[i];
        p = s;
    }
    if (p < 0) p = 0; if (p > 1) p = 1;
    Safefree(sc); ft_free(&S);
    return p;
}

static void calculate_exact_stats(long a, long b, long c, long d, double conf,
                                  const char *alt, double *orp, double *lop, double *hip) {
    ft_support S;
    if (!ft_init(&S, a, b, c, d)) { *orp = NAN; *lop = NAN; *hip = NAN; return; }
    double *sc; Newx(sc, S.ns, double);
    long x = S.x, lo = S.lo, hi = S.hi;

    /* conditional MLE of the odds ratio */
    double est;
    if      (x == lo) est = 0.0;
    else if (x == hi) est = INFINITY;
    else {
        double mu = ft_mnhyper(&S, 1.0, sc);
        ft_rc r = { &S, (double)x, sc, 0 };
        if      (mu > x) { r.mode = 0; est = ft_zeroin(0, 1, ft_rootf, &r, FT_TOL, 1000); }
        else if (mu < x) { r.mode = 1; est = 1.0 / ft_zeroin(FT_EPS, 1, ft_rootf, &r, FT_TOL, 1000); }
        else             est = 1.0;
    }
    *orp = est;

    /* confidence interval via inversion of the noncentral hypergeometric */
    double clo, chi;
    ft_rc r = { &S, 0, sc, 0 };
    #define FT_NCP_L(alpha, dst) do {                                                    \
        if (x == lo) { dst = 0.0; } else {                                               \
            double p = ft_pnhyper(&S, x, 1.0, 1, sc);                                     \
            if (p > (alpha))      { r.mode = 4; r.target = (alpha); dst = ft_zeroin(0, 1, ft_rootf, &r, FT_TOL, 1000); } \
            else if (p < (alpha)) { r.mode = 5; r.target = (alpha); dst = 1.0 / ft_zeroin(FT_EPS, 1, ft_rootf, &r, FT_TOL, 1000); } \
            else dst = 1.0; } } while (0)
    #define FT_NCP_U(alpha, dst) do {                                                    \
        if (x == hi) { dst = INFINITY; } else {                                          \
            double p = ft_pnhyper(&S, x, 1.0, 0, sc);                                     \
            if (p < (alpha))      { r.mode = 2; r.target = (alpha); dst = ft_zeroin(0, 1, ft_rootf, &r, FT_TOL, 1000); } \
            else if (p > (alpha)) { r.mode = 3; r.target = (alpha); dst = 1.0 / ft_zeroin(FT_EPS, 1, ft_rootf, &r, FT_TOL, 1000); } \
            else dst = 1.0; } } while (0)

    if      (!strcmp(alt, "less"))    { clo = 0.0;            FT_NCP_U(1 - conf, chi); }
    else if (!strcmp(alt, "greater")) { FT_NCP_L(1 - conf, clo); chi = INFINITY; }
    else { double al = (1 - conf) / 2; FT_NCP_L(al, clo); FT_NCP_U(al, chi); }

    *lop = clo; *hip = chi;
    Safefree(sc); ft_free(&S);
}

/* small helper: fetch a nonnegative integer cell from an SV, with validation */
static long ft_cell(pTHX_ SV *sv, const char *what) {
    if (!sv || !SvOK(sv)) croak("fisher_test: %s is undef", what);
    if (!looks_like_number(sv)) croak("fisher_test: %s is not a number", what);
    IV v = SvIV(sv);
    if (v < 0) croak("fisher_test: %s must be nonnegative (got %" IVdf ")", what, v);
    return (long)v;
}

/* ----------------------------------------------------------------- XS body */
/* Place this inside your MODULE = Stats::LikeR  PACKAGE = Stats::LikeR block. */

SV* fisher_test(...)
CODE:
{
    if (items < 1) croak("fisher_test requires at least a data reference");

    SV *data_ref = ST(0);
    NV conf_level = 0.95;
    const char *alternative = "two.sided";

    for (int i = 1; i < items; i += 2) {
        if (i + 1 >= items) croak("fisher_test: odd number of named arguments");
        const char *key = SvPV_nolen(ST(i));
        SV *val = ST(i + 1);
        if (strEQ(key, "conf_level") || strEQ(key, "conf.level")) {
            conf_level = SvNV(val);
            if (!(conf_level > 0 && conf_level < 1))
                croak("fisher_test: conf_level must be between 0 and 1");
        } else if (strEQ(key, "alternative")) {
            alternative = SvPV_nolen(val);
            if (strNE(alternative, "two.sided") && strNE(alternative, "less") &&
                strNE(alternative, "greater"))
                croak("fisher_test: alternative must be 'two.sided', 'less' or 'greater'");
        } else {
            croak("fisher_test: unknown argument '%s'", key);
        }
    }

    if (!SvROK(data_ref)) croak("fisher_test requires a reference to a 2x2 Array or Hash");
    SV *deref = SvRV(data_ref);
    long a = 0, b = 0, c = 0, d = 0;

    if (SvTYPE(deref) == SVt_PVAV) {
        AV *outer = (AV *)deref;
        if (av_len(outer) != 1) croak("Outer array must have exactly 2 rows");
        SV **r1p = av_fetch(outer, 0, 0);
        SV **r2p = av_fetch(outer, 1, 0);
        if (!(r1p && r2p && SvROK(*r1p) && SvROK(*r2p)
              && SvTYPE(SvRV(*r1p)) == SVt_PVAV && SvTYPE(SvRV(*r2p)) == SVt_PVAV))
            croak("Invalid 2D array structure: need two array-ref rows");
        AV *r1 = (AV *)SvRV(*r1p), *r2 = (AV *)SvRV(*r2p);
        if (av_len(r1) != 1 || av_len(r2) != 1)
            croak("Each row must have exactly 2 columns");
        a = ft_cell(aTHX_ *av_fetch(r1, 0, 0), "cell [0][0]");
        b = ft_cell(aTHX_ *av_fetch(r1, 1, 0), "cell [0][1]");
        c = ft_cell(aTHX_ *av_fetch(r2, 0, 0), "cell [1][0]");
        d = ft_cell(aTHX_ *av_fetch(r2, 1, 0), "cell [1][1]");
    }
    else if (SvTYPE(deref) == SVt_PVHV) {
        /* 2x2 hash; rows and columns are ordered by lexical key sort so the
         * result is deterministic regardless of Perl's hash randomization. */
        HV *outer = (HV *)deref;
        if (HvUSEDKEYS(outer) != 2) croak("Outer hash must have exactly 2 keys");
        hv_iterinit(outer);
        HE *e1 = hv_iternext(outer), *e2 = hv_iternext(outer);
        const char *ok1 = SvPV_nolen(hv_iterkeysv(e1));
        int swap_rows = strcmp(ok1, SvPV_nolen(hv_iterkeysv(e2))) > 0;
        SV *row1_sv = hv_iterval(outer, swap_rows ? e2 : e1);
        SV *row2_sv = hv_iterval(outer, swap_rows ? e1 : e2);
        if (!SvROK(row1_sv) || SvTYPE(SvRV(row1_sv)) != SVt_PVHV ||
            !SvROK(row2_sv) || SvTYPE(SvRV(row2_sv)) != SVt_PVHV)
            croak("Inner elements must be hash refs");

        HV *rows[2]; rows[0] = (HV *)SvRV(row1_sv); rows[1] = (HV *)SvRV(row2_sv);
        long cells[2][2];
        for (int rr = 0; rr < 2; rr++) {
            HV *in = rows[rr];
            if (HvUSEDKEYS(in) != 2) croak("Inner hashes must have exactly 2 keys");
            hv_iterinit(in);
            HE *c1 = hv_iternext(in), *c2 = hv_iternext(in);
            const char *k1 = SvPV_nolen(hv_iterkeysv(c1));
            int swap_cols = strcmp(k1, SvPV_nolen(hv_iterkeysv(c2))) > 0;
            HE *col0 = swap_cols ? c2 : c1;
            HE *col1 = swap_cols ? c1 : c2;
            cells[rr][0] = ft_cell(aTHX_ hv_iterval(in, col0), "hash cell");
            cells[rr][1] = ft_cell(aTHX_ hv_iterval(in, col1), "hash cell");
        }
        a = cells[0][0]; b = cells[0][1]; c = cells[1][0]; d = cells[1][1];
    }
    else {
        croak("Input must be a 2D Array or 2D Hash");
    }

    if (a + b + c + d == 0) croak("fisher_test: table is all zeros");

    NV p_val = exact_p_value(a, b, c, d, alternative);
    NV mle_or, ci_low, ci_high;
    calculate_exact_stats(a, b, c, d, conf_level, alternative, &mle_or, &ci_low, &ci_high);

    HV *ret = newHV();
    hv_stores(ret, "method", newSVpv("Fisher's Exact Test for Count Data", 0));
    hv_stores(ret, "alternative", newSVpv(alternative, 0));
    AV *ci = newAV();
    av_push(ci, newSVnv(ci_low));
    av_push(ci, newSVnv(ci_high));
    hv_stores(ret, "conf_int", newRV_noinc((SV *)ci));
    HV *est = newHV();
    hv_stores(est, "odds ratio", newSVnv(mle_or));
    hv_stores(ret, "estimate", newRV_noinc((SV *)est));
    hv_stores(ret, "p_value", newSVnv(p_val));
    hv_stores(ret, "conf_level", newSVnv(conf_level));
    RETVAL = newRV_noinc((SV *)ret);
}
OUTPUT:
  RETVAL
