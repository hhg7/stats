#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"
#include <math.h>
#include <string.h>

/* --- C HELPER SECTION --- */
#include <ctype.h>
#include <stdlib.h>

/* Struct for sorting p-values while remembering their original index */
typedef struct {
    double p;
    unsigned int orig_idx;
} PVal;

/* Comparator for qsort */
static int cmp_pval(const void *a, const void *b) {
    double diff = ((PVal*)a)->p - ((PVal*)b)->p;
    if (diff < 0) return -1;
    if (diff > 0) return 1;
    /* Stabilize sort by falling back to original index */
    return ((PVal*)a)->orig_idx - ((PVal*)b)->orig_idx; 
}

/* Comparator for sorting raw doubles (used by median) */
static int cmp_double(const void *a, const void *b) {
    double da = *(const double *)a;
    double db = *(const double *)b;
    if (da < db) return -1;
    if (da > db) return  1;
    return 0;
}

/* -----------------------------------------------------------------------
 * Helpers for cor(): ranking (Spearman), Pearson r, Kendall tau-b
 * ----------------------------------------------------------------------- */

/* Item used to sort values while remembering their original index,
 * needed for average-rank tie-breaking in Spearman correlation.        */
typedef struct {
    double       val;
    unsigned int idx;
} RankItem;

static int cmp_rank_item(const void *a, const void *b) {
    double diff = ((RankItem*)a)->val - ((RankItem*)b)->val;
    if (diff < 0) return -1;
    if (diff > 0) return  1;
    return 0;
}

/* Compute 1-based average ranks with tie-breaking into out[].
 * in[] is not modified.                                                 */
static void rank_data(const double *in, double *out, unsigned int n) {
    RankItem *ri;
    Newx(ri, n, RankItem);
    for (unsigned int i = 0; i < n; i++) { ri[i].val = in[i]; ri[i].idx = i; }
    qsort(ri, n, sizeof(RankItem), cmp_rank_item);

    unsigned int i = 0;
    while (i < n) {
        unsigned int j = i;
        /* Find the full extent of this tie group */
        while (j + 1 < n && ri[j + 1].val == ri[j].val) j++;
        /* All members get the average of ranks i+1 … j+1 (1-based) */
        double avg = (double)(i + j) / 2.0 + 1.0;
        for (unsigned int k = i; k <= j; k++) out[ri[k].idx] = avg;
        i = j + 1;
    }
    Safefree(ri);
}

/* Pearson product-moment r between two n-element arrays.
 * Returns NAN when either variable has zero variance (matches R).       */
static double pearson_corr(const double *x, const double *y, unsigned int n) {
    double sx = 0, sy = 0, sxy = 0, sx2 = 0, sy2 = 0;
    for (unsigned int i = 0; i < n; i++) {
        sx  += x[i];     sy  += y[i];
        sxy += x[i]*y[i]; sx2 += x[i]*x[i]; sy2 += y[i]*y[i];
    }
    double num = (double)n * sxy - sx * sy;
    double den = sqrt(((double)n * sx2 - sx*sx) * ((double)n * sy2 - sy*sy));
    if (den == 0.0) return NAN;
    return num / den;
}

/* Kendall's tau-b between two n-element arrays.
 *
 *   tau-b = (C − D) / sqrt((C + D + T_x)(C + D + T_y))
 *
 * where C = concordant pairs, D = discordant, T_x = pairs tied only on
 * x, T_y = pairs tied only on y.  Joint ties (both zero) are excluded
 * from numerator and denominator, matching R's cor(method="kendall").
 * Returns NAN when the denominator is zero.                             */
static double kendall_tau_b(const double *x, const double *y, unsigned int n) {
    long long C = 0, D = 0, tie_x = 0, tie_y = 0;
    for (unsigned int i = 0; i < n - 1; i++) {
        for (unsigned int j = i + 1; j < n; j++) {
            int sx = (x[i] > x[j]) - (x[i] < x[j]);   /* sign of x[i]-x[j] */
            int sy = (y[i] > y[j]) - (y[i] < y[j]);
            if      (sx == 0 && sy == 0) { /* joint tie — not counted */ }
            else if (sx == 0)            tie_x++;
            else if (sy == 0)            tie_y++;
            else if (sx == sy)           C++;
            else                         D++;
        }
    }
    double denom = sqrt((double)(C + D + tie_x) * (double)(C + D + tie_y));
    if (denom == 0.0) return NAN;
    return (double)(C - D) / denom;
}

/* Single dispatch: compute correlation according to method string.
 * Allocates and frees temporary rank arrays internally for Spearman.   */
static double compute_cor(const double *x, const double *y,
                           unsigned int n, const char *method) {
    if (strcmp(method, "spearman") == 0) {
        double *rx, *ry;
        Newx(rx, n, double); Newx(ry, n, double);
        rank_data(x, rx, n);
        rank_data(y, ry, n);
        double r = pearson_corr(rx, ry, n);
        Safefree(rx); Safefree(ry);
        return r;
    }
    if (strcmp(method, "kendall") == 0)
        return kendall_tau_b(x, y, n);
    /* default: pearson */
    return pearson_corr(x, y, n);
}

/* Math macros */
#define MAX_ITER 200
#define EPS 3.0e-7
#define FPMIN 1.0e-30

static double _incbeta_cf(double a, double b, double x) {
    int m;
    double aa, c, d, del, h, qab, qam, qap;
    qab = a + b; qap = a + 1.0; qam = a - 1.0;
    c = 1.0; d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d; h = d;
    for (m = 1; m <= MAX_ITER; m++) {
        int m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d; del = d * c; h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    return h;
}

static double incbeta(double a, double b, double x) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    double bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0)) return bt * _incbeta_cf(a, b, x) / a;
    return 1.0 - bt * _incbeta_cf(b, a, 1.0 - x) / b;
}

static double get_t_pvalue(double t, double df, const char* alt) {
    double x = df / (df + t * t);
    double prob_2tail = incbeta(df / 2.0, 0.5, x);
    if (strcmp(alt, "less") == 0) return (t < 0) ? 0.5 * prob_2tail : 1.0 - 0.5 * prob_2tail;
    if (strcmp(alt, "greater") == 0) return (t > 0) ? 0.5 * prob_2tail : 1.0 - 0.5 * prob_2tail;
    return prob_2tail;
}

/* Bisection algorithm to find the inverse t-distribution (Critical t-value) */
static double qt_tail(double df, double p_tail) {
    double low = 0.0, high = 1.0;
    /* Find upper bound */
    while (get_t_pvalue(high, df, "greater") > p_tail) {
        low = high;
        high *= 2.0;
        if (high > 1000000.0) break; /* Fallback limit */
    }
    /* Bisect to find the root */
    for (unsigned int i = 0; i < 100; i++) {
        double mid = (low + high) / 2.0;
        double p_mid = get_t_pvalue(mid, df, "greater");
        if (p_mid > p_tail) {
            low = mid;
        } else {
            high = mid;
        }
        if (high - low < 1e-8) break;
    }
    return (low + high) / 2.0;
}

/* --- XS SECTION --- */

MODULE = stats		PACKAGE = stats		

double
mean(...)
    PROTOTYPE: @
    INIT:
        double total = 0;
        unsigned int count = 0;
    CODE:
        for (unsigned int i = 0; i < items; i++) {
            SV* arg = ST(i);
            if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
                AV* av = (AV*)SvRV(arg);
                unsigned int len = av_len(av) + 1;
                for (unsigned int j = 0; j < len; j++) {
                    SV** tv = av_fetch(av, j, 0);
                    if (tv && SvOK(*tv)) { total += SvNV(*tv); count++; }
                }
            } else if (SvOK(arg)) {
                total += SvNV(arg); count++;
            }
        }
        if (count == 0) croak("mean needs >= 1 element");
        RETVAL = total / count;
    OUTPUT:
        RETVAL

double
stdev(...)
    PROTOTYPE: @
    INIT:
        double total = 0;
        unsigned int count = 0;
        double sum_sq = 0;
    CODE:
        for (unsigned int i = 0; i < items; i++) {
            SV* arg = ST(i);
            if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
                AV* av = (AV*)SvRV(arg);
                unsigned int len = av_len(av) + 1;
                for (unsigned int j = 0; j < len; j++) {
                    SV** tv = av_fetch(av, j, 0);
                    if (tv && SvOK(*tv)) { total += SvNV(*tv); count++; }
                }
            } else if (SvOK(arg)) {
                total += SvNV(arg); count++;
            }
        }
        if (count < 2) croak("stdev needs >= 2 elements");
        double avg = total / count;
        for (unsigned int i = 0; i < items; i++) {
            SV* arg = ST(i);
            if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
                AV* av = (AV*)SvRV(arg);
                unsigned int len = av_len(av) + 1;
                for (unsigned int j = 0; j < len; j++) {
                    SV** tv = av_fetch(av, j, 0);
                    if (tv && SvOK(*tv)) sum_sq += pow(SvNV(*tv) - avg, 2);
                }
            } else if (SvOK(arg)) {
                sum_sq += pow(SvNV(arg) - avg, 2);
            }
        }
        RETVAL = sqrt(sum_sq / (count - 1));
    OUTPUT:
        RETVAL

double
pearson_r(SV* x_sv, SV* y_sv)
    INIT:
        if (!SvROK(x_sv) || SvTYPE(SvRV(x_sv)) != SVt_PVAV ||
            !SvROK(y_sv) || SvTYPE(SvRV(y_sv)) != SVt_PVAV) {
            croak("Both arguments must be array references");
        }
        AV* x_av = (AV*)SvRV(x_sv);
        AV* y_av = (AV*)SvRV(y_sv);
        unsigned int n = av_len(x_av) + 1;
        if (n != (av_len(y_av) + 1)) croak("Arrays must have the same number of elements");
        if (n < 2) croak("Need at least 2 elements");
        double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
    CODE:
        for (unsigned int i = 0; i < n; i++) {
            SV** xv = av_fetch(x_av, i, 0);
            SV** yv = av_fetch(y_av, i, 0);
            double x = (xv && SvOK(*xv)) ? SvNV(*xv) : 0.0;
            double y = (yv && SvOK(*yv)) ? SvNV(*yv) : 0.0;
            sum_x += x; sum_y += y; sum_xy += x * y;
            sum_x2 += x * x; sum_y2 += y * y;
        }
        double num = (n * sum_xy) - (sum_x * sum_y);
        double den = sqrt((n * sum_x2 - (sum_x * sum_x)) * (n * sum_y2 - (sum_y * sum_y)));
        if (den == 0) XSRETURN_UNDEF;
        RETVAL = num / den;
    OUTPUT:
        RETVAL


SV* t_test(SV* args_sv)
    INIT:
        if (!SvROK(args_sv) || SvTYPE(SvRV(args_sv)) != SVt_PVHV) {
            croak("Usage: t_test({ x => [...], y => [...], ... }) - Arguments must be a HASH reference");
        }
        HV* args_hv = (HV*)SvRV(args_sv);
        SV** temp;

        if (!(temp = hv_fetch(args_hv, "x", 1, 0)) || !SvROK(*temp) || SvTYPE(SvRV(*temp)) != SVt_PVAV) {
            croak("t_test: 'x' is a required argument and must be an ARRAY reference");
        }
        AV* x_av = (AV*)SvRV(*temp);
        unsigned int nx = av_len(x_av) + 1;
        if (nx < 2) croak("t_test: 'x' needs at least 2 elements");

        AV* y_av = NULL;
        if ((temp = hv_fetch(args_hv, "y", 1, 0)) && SvROK(*temp) && SvTYPE(SvRV(*temp)) == SVt_PVAV) {
            y_av = (AV*)SvRV(*temp);
        }

        double mu = 0.0;
        if ((temp = hv_fetch(args_hv, "mu", 2, 0))) mu = SvNV(*temp);

        bool paired = FALSE;
        if ((temp = hv_fetch(args_hv, "paired", 6, 0))) paired = SvTRUE(*temp);

        bool var_equal = FALSE;
        if ((temp = hv_fetch(args_hv, "var_equal", 9, 0))) var_equal = SvTRUE(*temp);

        double conf_level = 0.95;
        if ((temp = hv_fetch(args_hv, "conf_level", 10, 0))) conf_level = SvNV(*temp);
        if (conf_level <= 0.0 || conf_level >= 1.0) croak("t_test: 'conf_level' must be between 0 and 1");

        const char* alternative = "two.sided";
        if ((temp = hv_fetch(args_hv, "alternative", 11, 0))) alternative = SvPV_nolen(*temp);

        double sum_x = 0, sum_x2 = 0, mean_x, var_x;
        double t_stat, df, p_val, std_err, cint_est;
        HV* results = newHV();

    CODE:
        for (unsigned int i = 0; i < nx; i++) {
            SV** tv = av_fetch(x_av, i, 0);
            double val = (tv && SvOK(*tv)) ? SvNV(*tv) : 0;
            sum_x += val; sum_x2 += val * val;
        }
        mean_x = sum_x / nx;
        var_x = (sum_x2 - (sum_x * sum_x) / nx) / (nx - 1);

        if (paired || y_av) {
            if (!y_av) croak("t_test: 'y' must be provided for paired or two-sample tests");
            unsigned int ny = av_len(y_av) + 1;
            if (paired && ny != nx) croak("t_test: Paired arrays must be same length");

            double sum_y = 0, sum_y2 = 0, mean_y, var_y;
            for (unsigned int i = 0; i < ny; i++) {
                SV** tv = av_fetch(y_av, i, 0);
                double val = (tv && SvOK(*tv)) ? SvNV(*tv) : 0;
                sum_y += val; sum_y2 += val * val;
            }
            mean_y = sum_y / ny;
            var_y = (sum_y2 - (sum_y * sum_y) / ny) / (ny - 1);

            if (paired) {
                double sum_d = 0, sum_d2 = 0;
                for (unsigned int i = 0; i < nx; i++) {
                    double dx = SvNV(*av_fetch(x_av, i, 0));
                    double dy = SvNV(*av_fetch(y_av, i, 0));
                    sum_d += (dx - dy); sum_d2 += (dx - dy) * (dx - dy);
                }
                double mean_d = sum_d / nx;
                double var_d = (sum_d2 - (sum_d * sum_d) / nx) / (nx - 1);
                
                cint_est = mean_d;
                std_err = sqrt(var_d / nx);
                t_stat = (cint_est - mu) / std_err;
                df = nx - 1;
                hv_store(results, "estimate", 8, newSVnv(mean_d), 0);
            } else if (var_equal) {
                double pooled_var = ((nx - 1) * var_x + (ny - 1) * var_y) / (nx + ny - 2);
                cint_est = mean_x - mean_y;
                std_err = sqrt(pooled_var * (1.0 / nx + 1.0 / ny));
                t_stat = (cint_est - mu) / std_err;
                df = nx + ny - 2;
                hv_store(results, "estimate_x", 10, newSVnv(mean_x), 0);
                hv_store(results, "estimate_y", 10, newSVnv(mean_y), 0);
            } else {
                cint_est = mean_x - mean_y;
                double stderr_x2 = var_x / nx;
                double stderr_y2 = var_y / ny;
                std_err = sqrt(stderr_x2 + stderr_y2);
                t_stat = (cint_est - mu) / std_err;
                df = pow(stderr_x2 + stderr_y2, 2) / (pow(stderr_x2, 2) / (nx - 1) + pow(stderr_y2, 2) / (ny - 1));
                hv_store(results, "estimate_x", 10, newSVnv(mean_x), 0);
                hv_store(results, "estimate_y", 10, newSVnv(mean_y), 0);
            }
        } else {
            cint_est = mean_x;
            std_err = sqrt(var_x / nx);
            t_stat = (cint_est - mu) / std_err;
            df = nx - 1;
            hv_store(results, "estimate", 8, newSVnv(mean_x), 0);
        }

        /* Calculate P-Value */
        p_val = get_t_pvalue(t_stat, df, alternative);
        
        /* Calculate Confidence Intervals */
        double alpha = 1.0 - conf_level;
        double t_crit = 0.0;
        double ci_lower, ci_upper;

        if (strcmp(alternative, "less") == 0) {
            t_crit = qt_tail(df, alpha);
            ci_lower = -INFINITY;
            ci_upper = cint_est + t_crit * std_err;
        } else if (strcmp(alternative, "greater") == 0) {
            t_crit = qt_tail(df, alpha);
            ci_lower = cint_est - t_crit * std_err;
            ci_upper = INFINITY;
        } else {
            t_crit = qt_tail(df, alpha / 2.0);
            ci_lower = cint_est - t_crit * std_err;
            ci_upper = cint_est + t_crit * std_err;
        }

        AV* conf_int = newAV();
        av_push(conf_int, newSVnv(ci_lower));
        av_push(conf_int, newSVnv(ci_upper));

        hv_store(results, "statistic", 9, newSVnv(t_stat), 0);
        hv_store(results, "df", 2, newSVnv(df), 0);
        hv_store(results, "p_value", 7, newSVnv(p_val), 0);
        hv_store(results, "conf_int", 8, newRV_noinc((SV*)conf_int), 0);
        
        RETVAL = newRV_noinc((SV*)results);
    OUTPUT:
        RETVAL

void
p_adjust(SV* p_sv, const char* method = "holm")
    INIT:
        if (!SvROK(p_sv) || SvTYPE(SvRV(p_sv)) != SVt_PVAV) {
            croak("p_adjust: first argument must be an ARRAY reference of p-values");
        }
        
        AV* p_av = (AV*)SvRV(p_sv);
        unsigned int n = av_len(p_av) + 1;
        
        /* Handle empty input */
        if (n == 0) {
            XSRETURN_EMPTY;
        }

        /* Normalize method string */
        char meth[64];
        strncpy(meth, method, 63); meth[63] = '\0';
        for(int i = 0; meth[i]; i++) meth[i] = tolower(meth[i]);

        /* Resolve aliases */
        if (strstr(meth, "benjamini") && strstr(meth, "hochberg")) strcpy(meth, "bh");
        if (strstr(meth, "benjamini") && strstr(meth, "yekutieli")) strcpy(meth, "by");
        if (strcmp(meth, "fdr") == 0) strcpy(meth, "bh");

        /* Allocate C memory */
        PVal* arr;
        double* adj;
        Newx(arr, n, PVal);
        Newx(adj, n, double);

        for (unsigned int i = 0; i < n; i++) {
            SV** tv = av_fetch(p_av, i, 0);
            arr[i].p = (tv && SvOK(*tv)) ? SvNV(*tv) : 1.0;
            arr[i].orig_idx = i;
        }

        /* Sort ascending (Stable sort using original index) */
        qsort(arr, n, sizeof(PVal), cmp_pval);

    PPCODE:
        if (strcmp(meth, "bonferroni") == 0) {
            for (int i = 0; i < n; i++) {
                double v = arr[i].p * n;
                adj[arr[i].orig_idx] = (v < 1.0) ? v : 1.0;
            }
        } else if (strcmp(meth, "holm") == 0) {
            double cummax = 0.0;
            for (int i = 0; i < n; i++) {
                double v = arr[i].p * (n - i);
                if (v > cummax) cummax = v;
                adj[arr[i].orig_idx] = (cummax < 1.0) ? cummax : 1.0;
            }
        } else if (strcmp(meth, "hochberg") == 0) {
            double cummin = 1.0;
            for (int i = n - 1; i >= 0; i--) {
                double v = arr[i].p * (n - i);
                if (v < cummin) cummin = v;
                adj[arr[i].orig_idx] = (cummin < 1.0) ? cummin : 1.0;
            }
        } else if (strcmp(meth, "bh") == 0) {
            double cummin = 1.0;
            for (int i = n - 1; i >= 0; i--) {
                double v = arr[i].p * n / (i + 1.0);
                if (v < cummin) cummin = v;
                adj[arr[i].orig_idx] = (cummin < 1.0) ? cummin : 1.0;
            }
        } else if (strcmp(meth, "by") == 0) {
            double q = 0.0;
            for (int i = 1; i <= n; i++) q += 1.0 / i;
            double cummin = 1.0;
            for (int i = n - 1; i >= 0; i--) {
                double v = arr[i].p * n / (i + 1.0) * q;
                if (v < cummin) cummin = v;
                adj[arr[i].orig_idx] = (cummin < 1.0) ? cummin : 1.0;
            }
        } else if (strcmp(meth, "hommel") == 0) {
            double *pa = (double *)malloc(n * sizeof(double));
            double *q_arr = (double *)malloc(n * sizeof(double));
            
            /* Initial: min(n * p[i] / (i + 1)) */
            double min_val = n * arr[0].p;
            for (unsigned int i = 1; i < n; i++) {
                double temp = (n * arr[i].p) / (i + 1.0);
                if (temp < min_val) {
                    min_val = temp;
                }
            }
            /* pa <- q <- rep(min, n) */
            for (unsigned int i = 0; i < n; i++) {
                pa[i] = min_val;
                q_arr[i] = min_val;
            }
            for (unsigned int j = n - 1; j >= 2; j--) {
                int n_mj = n - j;       /* Max index for 'ij'. Length is n_mj + 1 */
                int i2_len = j - 1;     /* Length of 'i2' */
                
                /* Calculate q1 = min(j * p[i2] / (2:j)) */
                double q1 = (j * arr[n_mj + 1].p) / 2.0;
                for (unsigned int k = 1; k < i2_len; k++) {
                    double temp_q1 = (j * arr[n_mj + 1 + k].p) / (2.0 + k);
                    if (temp_q1 < q1) {
                        q1 = temp_q1;
                    }
                }

                /* q[ij] <- pmin(j * p[ij], q1) */
                for (unsigned int i = 0; i <= n_mj; i++) {
                    double v = j * arr[i].p;
                    q_arr[i] = (v < q1) ? v : q1;
                }

                /* q[i2] <- q[n - j] */
                for (unsigned int i = 0; i < i2_len; i++) {
                    q_arr[n_mj + 1 + i] = q_arr[n_mj];
                }

                /* pa <- pmax(pa, q) */
                for (unsigned int i = 0; i < n; i++) {
                    if (pa[i] < q_arr[i]) {
                        pa[i] = q_arr[i];
                    }
                }
            }
            
            /* pmin(1, pmax(pa, p))[ro] — map sorted results back to original indices */
            for (unsigned int i = 0; i < n; i++) {
                double v = (pa[i] > arr[i].p) ? pa[i] : arr[i].p;
                if (v > 1.0) v = 1.0;
                adj[arr[i].orig_idx] = v;
            }
            free(pa); pa = NULL;
            free(q_arr); q_arr = NULL;
        } else if (strcmp(meth, "none") == 0) {
            for (unsigned int i = 0; i < n; i++) {
            	adj[arr[i].orig_idx] = arr[i].p;
            }
        } else {
            Safefree(arr); Safefree(adj);
            croak("Unknown p-value adjustment method: %s", method);
        }
        /* Push values onto the Perl stack as a flat list */
        EXTEND(SP, n);
        for (int i = 0; i < n; i++) {
            PUSHs(sv_2mortal(newSVnv(adj[i])));
        }

        Safefree(arr);
        Safefree(adj);


/* -----------------------------------------------------------------------
 * median(LIST | \@array ...)
 *
 * Matches base R median():
 *   - Accepts the same flat-list / array-ref input style as mean().
 *   - Returns the middle value for odd n, or the arithmetic mean of the
 *     two middle values for even n (R's default "type-7" behaviour).
 *   - Croaks on empty input.
 * ----------------------------------------------------------------------- */
double
median(...)
    PROTOTYPE: @
    INIT:
        unsigned int count = 0;
        unsigned int alloc = 64;
        double *vals;
        Newx(vals, alloc, double);
    CODE:
        /* Collect every numeric value from the argument list */
        for (unsigned int i = 0; i < items; i++) {
            SV* arg = ST(i);
            if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
                AV* av = (AV*)SvRV(arg);
                unsigned int len = av_len(av) + 1;
                for (unsigned int j = 0; j < len; j++) {
                    SV** tv = av_fetch(av, j, 0);
                    if (tv && SvOK(*tv)) {
                        if (count >= alloc) { alloc *= 2; Renew(vals, alloc, double); }
                        vals[count++] = SvNV(*tv);
                    }
                }
            } else if (SvOK(arg)) {
                if (count >= alloc) { alloc *= 2; Renew(vals, alloc, double); }
                vals[count++] = SvNV(arg);
            }
        }
        if (count == 0) { Safefree(vals); croak("median needs >= 1 element"); }

        qsort(vals, count, sizeof(double), cmp_double);

        if (count % 2 == 1)
            RETVAL = vals[count / 2];
        else
            RETVAL = (vals[count / 2 - 1] + vals[count / 2]) / 2.0;

        Safefree(vals);
    OUTPUT:
        RETVAL

/* -----------------------------------------------------------------------
 * var(LIST | \@array ...)
 *
 * Matches base R var():
 *   - Same flexible input style as mean() / stdev().
 *   - Returns the *sample* variance (denominator n-1), identical to R.
 *   - Croaks when fewer than 2 elements are supplied.
 * ----------------------------------------------------------------------- */
double
var(...)
    PROTOTYPE: @
    INIT:
        double total = 0;
        unsigned int count = 0;
        double sum_sq = 0;
    CODE:
        /* First pass: accumulate sum for the mean */
        for (unsigned int i = 0; i < items; i++) {
            SV* arg = ST(i);
            if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
                AV* av = (AV*)SvRV(arg);
                unsigned int len = av_len(av) + 1;
                for (unsigned int j = 0; j < len; j++) {
                    SV** tv = av_fetch(av, j, 0);
                    if (tv && SvOK(*tv)) { total += SvNV(*tv); count++; }
                }
            } else if (SvOK(arg)) {
                total += SvNV(arg); count++;
            }
        }
        if (count < 2) croak("var needs >= 2 elements");
        double avg = total / count;

        /* Second pass: sum of squared deviations */
        for (unsigned int i = 0; i < items; i++) {
            SV* arg = ST(i);
            if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
                AV* av = (AV*)SvRV(arg);
                unsigned int len = av_len(av) + 1;
                for (unsigned int j = 0; j < len; j++) {
                    SV** tv = av_fetch(av, j, 0);
                    if (tv && SvOK(*tv)) sum_sq += pow(SvNV(*tv) - avg, 2);
                }
            } else if (SvOK(arg)) {
                sum_sq += pow(SvNV(arg) - avg, 2);
            }
        }
        RETVAL = sum_sq / (count - 1);
    OUTPUT:
        RETVAL

/* -----------------------------------------------------------------------
 * scale(\@x, center = 1, scale = 1)
 *
 * Matches base R scale(x, center = TRUE, scale = TRUE):
 *   x      : array reference of numeric values (required).
 *   center : boolean (default 1/TRUE).  When true, subtract the column
 *            mean so the result has mean ~= 0.
 *   scale  : boolean (default 1/TRUE).  When true, divide by the sample
 *            SD so the result has SD ~= 1.  A zero-SD column is left
 *            unchanged, matching R's behaviour.
 *
 * Returns a hash reference:
 *   {
 *     x      => [ ...scaled values... ],
 *     center => $mean_used,   # undef when center was FALSE
 *     scale  => $sd_used,     # undef when scale was FALSE or SD == 0
 *   }
 * This mirrors the "scaled:center" / "scaled:scale" attributes that R
 * attaches to the matrix it returns.
 * ----------------------------------------------------------------------- */
SV*
scale(SV* x_sv, bool do_center = 1, bool do_scale = 1)
    INIT:
        if (!SvROK(x_sv) || SvTYPE(SvRV(x_sv)) != SVt_PVAV)
            croak("scale: first argument must be an ARRAY reference");
        AV* x_av = (AV*)SvRV(x_sv);
        unsigned int n = av_len(x_av) + 1;
        if (n == 0) croak("scale: input array must be non-empty");
        if (do_scale && n < 2)
            croak("scale: cannot compute SD with fewer than 2 elements");
    CODE:
        double *vals;
        Newx(vals, n, double);
        for (unsigned int i = 0; i < n; i++) {
            SV** tv = av_fetch(x_av, i, 0);
            vals[i] = (tv && SvOK(*tv)) ? SvNV(*tv) : 0.0;
        }

        /* --- Center: subtract the mean --- */
        double center_val = 0.0;
        if (do_center) {
            double sum = 0.0;
            for (unsigned int i = 0; i < n; i++) sum += vals[i];
            center_val = sum / n;
            for (unsigned int i = 0; i < n; i++) vals[i] -= center_val;
        }

        /* --- Scale: divide by sample SD ---
         * R computes the SD from the *centred* values when center is also
         * TRUE, so vals[] is already centred here; mean2 stays 0.       */
        double scale_val = 0.0;
        if (do_scale) {
            double mean2 = 0.0;
            if (!do_center) {
                double sum = 0.0;
                for (unsigned int i = 0; i < n; i++) sum += vals[i];
                mean2 = sum / n;
            }
            double ss = 0.0;
            for (unsigned int i = 0; i < n; i++) ss += pow(vals[i] - mean2, 2);
            double sd = sqrt(ss / (n - 1));

            /* Zero-SD column: leave unchanged, return undef for scale attr */
            if (sd > 0.0) {
                scale_val = sd;
                for (unsigned int i = 0; i < n; i++) vals[i] /= sd;
            }
        }

        /* --- Build return hashref --- */
        HV* ret = newHV();

        AV* out_av = newAV();
        for (unsigned int i = 0; i < n; i++)
            av_push(out_av, newSVnv(vals[i]));
        hv_store(ret, "x",      1, newRV_noinc((SV*)out_av), 0);

        /* center attr: the mean subtracted, or undef */
        hv_store(ret, "center", 6, do_center              ? newSVnv(center_val) : newSV(0), 0);

        /* scale attr: the SD divided by, or undef when disabled / zero-SD */
        hv_store(ret, "scale",  5, (do_scale && scale_val > 0.0) ? newSVnv(scale_val) : newSV(0), 0);

        Safefree(vals);
        RETVAL = newRV_noinc((SV*)ret);
    OUTPUT:
        RETVAL

/* -----------------------------------------------------------------------
 * cor(x, y = undef, method = "pearson")
 *
 * Matches base R cor(x, y = NULL, use = "everything", method = "pearson"):
 *
 *   x      : array ref of numbers  — a single variable (vector)
 *           OR array ref of array refs — a data matrix (rows × cols)
 *   y      : same shapes as x, or undef/omitted.
 *   method : "pearson" (default), "spearman", or "kendall"
 *
 * Return value depends on inputs (mirroring R):
 *
 *   cor(\@vec)              → 1 (scalar SV)
 *   cor(\@vec, \@vec)       → scalar r / tau
 *   cor(\@matrix)           → p×p symmetric AoA ref   (diagonal = 1)
 *   cor(\@matrix_x, \@matrix_y) → p×q cross-corr AoA ref
 *
 * Matrix AoA convention: outer array = rows, inner arrays = columns,
 * matching standard Perl idiom.  Correlation is between *columns*,
 * matching R's treatment of variables as columns.
 *
 * Off-diagonal entries of a symmetric result exploit symmetry
 * (computed once, mirrored) to halve the work.
 * ----------------------------------------------------------------------- */
SV*
cor(SV* x_sv, SV* y_sv = &PL_sv_undef, const char* method = "pearson")
    INIT:
        /* --- validate method ------------------------------------------- */
        if (strcmp(method, "pearson")  != 0 &&
            strcmp(method, "spearman") != 0 &&
            strcmp(method, "kendall")  != 0)
            croak("cor: unknown method '%s' (use 'pearson', 'spearman', or 'kendall')",
                  method);

        /* --- validate x ------------------------------------------------ */
        if (!SvROK(x_sv) || SvTYPE(SvRV(x_sv)) != SVt_PVAV)
            croak("cor: x must be an ARRAY reference");
        AV* x_av = (AV*)SvRV(x_sv);
        int nx   = av_len(x_av) + 1;
        if (nx == 0) croak("cor: x is empty");

        /* --- detect whether x is a flat vector or a matrix (AoA) ------- */
        int x_is_matrix = 0;
        {
            SV** fp = av_fetch(x_av, 0, 0);
            if (fp && SvROK(*fp) && SvTYPE(SvRV(*fp)) == SVt_PVAV)
                x_is_matrix = 1;
        }

        /* --- detect y -------------------------------------------------- */
        int has_y = (SvOK(y_sv) && SvROK(y_sv) &&
                     SvTYPE(SvRV(y_sv)) == SVt_PVAV);
        AV* y_av     = has_y ? (AV*)SvRV(y_sv) : NULL;
        int ny       = has_y ? av_len(y_av) + 1 : 0;
        int y_is_matrix = 0;
        if (has_y && ny > 0) {
            SV** fp = av_fetch(y_av, 0, 0);
            if (fp && SvROK(*fp) && SvTYPE(SvRV(*fp)) == SVt_PVAV)
                y_is_matrix = 1;
        }
    CODE:
        /* ================================================================
         * Branch 1: both inputs are flat vectors  →  scalar result
         * ================================================================ */
        if (!x_is_matrix && !y_is_matrix) {
            if (!has_y) {
                /* cor(vector) == 1 by definition */
                RETVAL = newSVnv(1.0);
            } else {
                if (nx != ny)
                    croak("cor: x and y must have the same length (%d vs %d)",
                          nx, ny);
                if (nx < 2)
                    croak("cor: need at least 2 observations");

                double *xd, *yd;
                Newx(xd, nx, double);
                Newx(yd, ny, double);

                for (int i = 0; i < nx; i++) {
                    SV** tv = av_fetch(x_av, i, 0);
                    xd[i] = (tv && SvOK(*tv)) ? SvNV(*tv) : 0.0;
                }
                for (int i = 0; i < ny; i++) {
                    SV** tv = av_fetch(y_av, i, 0);
                    yd[i] = (tv && SvOK(*tv)) ? SvNV(*tv) : 0.0;
                }

                double r = compute_cor(xd, yd, nx, method);
                Safefree(xd); Safefree(yd);
                RETVAL = newSVnv(r);
            }
        }
        /* ================================================================
         * Branch 2: x is a matrix (or y is a matrix)  →  AoA result
         * ================================================================ */
        else {
            /* -- resolve x matrix dimensions ---------------------------- */
            if (!x_is_matrix)
                croak("cor: x must be a matrix (array ref of array refs) "
                      "when y is a matrix");

            SV** xr0 = av_fetch(x_av, 0, 0);
            if (!xr0 || !SvROK(*xr0) || SvTYPE(SvRV(*xr0)) != SVt_PVAV)
                croak("cor: each row of x must be an ARRAY reference");
            int ncols_x = av_len((AV*)SvRV(*xr0)) + 1;
            if (ncols_x == 0) croak("cor: x matrix has zero columns");
            int nrows   = nx;   /* observations */

            /* -- extract x columns -------------------------------------- */
            double **col_x;
            Newx(col_x, ncols_x, double*);
            for (int j = 0; j < ncols_x; j++) {
                Newx(col_x[j], nrows, double);
                for (int i = 0; i < nrows; i++) {
                    SV** rv = av_fetch(x_av, i, 0);
                    if (!rv || !SvROK(*rv))
                        croak("cor: x row %d is not an array ref", i);
                    AV*  row = (AV*)SvRV(*rv);
                    SV** cv  = av_fetch(row, j, 0);
                    col_x[j][i] = (cv && SvOK(*cv)) ? SvNV(*cv) : 0.0;
                }
            }

            /* -- resolve y: separate matrix or re-use x (symmetric) ---- */
            int    ncols_y;
            double **col_y   = NULL;
            int    symmetric = 0;   /* 1 = cor(X) — result is symmetric */

            if (has_y && y_is_matrix) {
                /* cross-correlation: X (nrows × p) vs Y (nrows × q)      */
                if (ny != nrows)
                    croak("cor: x and y must have the same number of rows "
                          "(%d vs %d)", nrows, ny);
                SV** yr0 = av_fetch(y_av, 0, 0);
                if (!yr0 || !SvROK(*yr0) || SvTYPE(SvRV(*yr0)) != SVt_PVAV)
                    croak("cor: each row of y must be an ARRAY reference");
                ncols_y = av_len((AV*)SvRV(*yr0)) + 1;
                if (ncols_y == 0) croak("cor: y matrix has zero columns");

                Newx(col_y, ncols_y, double*);
                for (int j = 0; j < ncols_y; j++) {
                    Newx(col_y[j], nrows, double);
                    for (int i = 0; i < nrows; i++) {
                        SV** rv = av_fetch(y_av, i, 0);
                        if (!rv || !SvROK(*rv))
                            croak("cor: y row %d is not an array ref", i);
                        AV*  row = (AV*)SvRV(*rv);
                        SV** cv  = av_fetch(row, j, 0);
                        col_y[j][i] = (cv && SvOK(*cv)) ? SvNV(*cv) : 0.0;
                    }
                }
            } else {
                /* cor(X) — symmetric p×p result; share column arrays     */
                ncols_y  = ncols_x;
                col_y    = col_x;
                symmetric = 1;
            }

            if (nrows < 2)
                croak("cor: need at least 2 observations (got %d)", nrows);

            /* -- build cache for symmetric case: compute upper triangle,
               store results, mirror to lower triangle                    */
            AV* result_av = newAV();
            av_extend(result_av, ncols_x - 1);

            /* Allocate per-row AVs up front so we can fill them in order */
            AV **rows_out;
            Newx(rows_out, ncols_x, AV*);
            for (int i = 0; i < ncols_x; i++) {
                rows_out[i] = newAV();
                av_extend(rows_out[i], ncols_y - 1);
            }

            if (symmetric) {
                /* Upper triangle + diagonal, then mirror.
                   r_cache[i][j] (j >= i) holds the computed value.      */
                double **r_cache;
                Newx(r_cache, ncols_x, double*);
                for (int i = 0; i < ncols_x; i++)
                    Newx(r_cache[i], ncols_x, double);

                for (int i = 0; i < ncols_x; i++) {
                    r_cache[i][i] = 1.0;                /* diagonal */
                    for (int j = i + 1; j < ncols_x; j++) {
                        double r = compute_cor(col_x[i], col_x[j], nrows, method);
                        r_cache[i][j] = r;
                        r_cache[j][i] = r;              /* symmetry */
                    }
                }
                /* fill output AoA from cache */
                for (int i = 0; i < ncols_x; i++)
                    for (int j = 0; j < ncols_x; j++)
                        av_store(rows_out[i], j, newSVnv(r_cache[i][j]));

                for (int i = 0; i < ncols_x; i++) Safefree(r_cache[i]);
                Safefree(r_cache);
            } else {
                /* cross-correlation: every (i,j) pair is independent     */
                for (int i = 0; i < ncols_x; i++)
                    for (int j = 0; j < ncols_y; j++)
                        av_store(rows_out[i], j,
                                 newSVnv(compute_cor(col_x[i], col_y[j],
                                                     nrows, method)));
            }

            /* push row AVs into result */
            for (int i = 0; i < ncols_x; i++)
                av_store(result_av, i, newRV_noinc((SV*)rows_out[i]));
            Safefree(rows_out);

            /* -- free column arrays ------------------------------------- */
            for (int j = 0; j < ncols_x; j++) Safefree(col_x[j]);
            Safefree(col_x);
            if (!symmetric) {
                for (int j = 0; j < ncols_y; j++) Safefree(col_y[j]);
                Safefree(col_y);
            }

            RETVAL = newRV_noinc((SV*)result_av);
        }
    OUTPUT:
        RETVAL
