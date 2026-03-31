/* --- C HELPER SECTION --- */
#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
/* Helper: log combination */
static double log_choose(size_t n, size_t k) {
    return lgamma((double)n + 1.0) - lgamma((double)k + 1.0) - lgamma((double)(n - k) + 1.0);
}


static double bisection_root(double (*func)(size_t, size_t, size_t, double),
                             size_t r1, size_t r2, size_t c1, double target,
                             double log_low0, double log_high0) {
    double log_low  = log_low0;
    double log_high = log_high0;
    double best_omega = exp((log_low + log_high) * 0.5);
    double best_error = 1e9;

    for (unsigned int i = 0; i < 20000; ++i) {
        double log_mid = 0.5 * (log_low + log_high);
        double omega   = exp(log_mid);
        double val     = func(r1, r2, c1, omega);
        double error   = fabs(val - target);

        if (error < best_error) {
            best_error = error;
            best_omega = omega;
        }

        if (val > target) {
            log_high = log_mid;
        } else {
            log_low  = log_mid;
        }
        if (log_high - log_low < 1e-30) break;
    }
    return best_omega;
}
/* Stable recursive-ratio PMF (matches R's internal method) */
static double dnoncentral_hypergeom(size_t x, size_t r1, size_t r2, size_t c1, double omega) {
    size_t min_x = (r2 > c1) ? 0 : c1 - r2;
    size_t max_x = (r1 < c1) ? r1 : c1;
    if (x < min_x || x > max_x || omega <= 0.0) return 0.0;

    double p = 1.0;      /* relative mass at min_x */
    double sum = 1.0;
    for (size_t i = min_x; i < x; ++i) {
        double ratio = ((double)(r1 - i) * (c1 - i) * omega) /
                       ((i + 1.0) * (r2 - c1 + i + 1.0));
        p *= ratio;
        sum += p;
    }
    return p / sum;
}

/* Expected value using the same stable recursive method */
static double expected_a(size_t r1, size_t r2, size_t c1, double omega) {
    size_t min_x = (r2 > c1) ? 0 : c1 - r2;
    size_t max_x = (r1 < c1) ? r1 : c1;
    if (omega <= 0.0) return (double)min_x;

    double p = 1.0;
    double sum_p = 1.0;
    double sum_xp = (double)min_x;
    for (size_t i = min_x; i < max_x; ++i) {
        double ratio = ((double)(r1 - i) * (c1 - i) * omega) /
                       ((i + 1.0) * (r2 - c1 + i + 1.0));
        p *= ratio;
        sum_p += p;
        sum_xp += (double)(i + 1) * p;
    }
    return sum_xp / sum_p;
}

/* Tail probabilities using recursive ratio (fast + stable, same as R) */
static void calc_tails(size_t a, size_t b, size_t c, size_t d, double omega, 
                       double *restrict lower_tail, double *restrict upper_tail) {
    size_t r1 = a + b, r2 = c + d, c1 = a + c;
    size_t min_x = (r2 > c1) ? 0 : c1 - r2;
    size_t max_x = (r1 < c1) ? r1 : c1;

    *lower_tail = 0.0;
    *upper_tail = 0.0;
    if (omega <= 0.0) {
        *lower_tail = 1.0;
        *upper_tail = 1.0;
        return;
    }

    double p = 1.0;
    double total = 1.0;
    if (min_x <= a) *lower_tail += p;
    if (min_x >= a) *upper_tail += p;

    for (size_t x = min_x; x < max_x; ++x) {
        double ratio = ((double)(r1 - x) * (c1 - x) * omega) /
                       ((x + 1.0) * (r2 - c1 + x + 1.0));
        p *= ratio;
        total += p;

        size_t curr = x + 1;
        if (curr <= a) *lower_tail += p;
        if (curr >= a) *upper_tail += p;
    }
    double inv = 1.0 / total;
    *lower_tail *= inv;
    *upper_tail *= inv;
}

/* Bisection with best-point tracking (matches R printed precision) */
static void calculate_exact_stats(size_t a, size_t b, size_t c, size_t d, double conf_level,
                                  double *restrict mle_or, double *restrict ci_low, double *restrict ci_high) {
    double alpha = 1.0 - conf_level;
    double target_alpha = alpha / 2.0;

    size_t r1 = a + b, r2 = c + d, c1 = a + c;
    size_t min_x = (r2 > c1) ? 0 : c1 - r2;
    size_t max_x = (r1 < c1) ? r1 : c1;
    /* MLE */
    if (a == min_x && a == max_x) *mle_or = 1.0;
    else if (a == min_x) *mle_or = 0.0;
    else if (a == max_x) *mle_or = INFINITY;
    else *mle_or = bisection_root(expected_a, r1, r2, c1, (double)a, -100.0, 100.0);
    /* Lower CI: P(X ≥ a | ω) = α/2 */
    if (a == min_x) {
        *ci_low = 0.0;
    } else {
        double log_low = -100.0, log_high = 100.0;
        double best = 1.0;
        double best_err = 1e9;
        double lt, ut;
        for (int i = 0; i < 10000; ++i) {
            double log_mid = 0.5 * (log_low + log_high);
            double mid = exp(log_mid);
            calc_tails(a, b, c, d, mid, &lt, &ut);
            double err = fabs(ut - target_alpha);
            if (err < best_err) { best_err = err; best = mid; }
            if (ut > target_alpha) log_high = log_mid; else log_low = log_mid;
            if (log_high - log_low < 1e-32) break;
        }
        *ci_low = best;
    }

    /* Upper CI: P(X ≤ a | ω) = α/2 */
    if (a == max_x) {
        *ci_high = INFINITY;
    } else {
        double log_low = -100.0, log_high = 100.0, best = 1.0, best_err = 1e9, lt, ut;
        for (int i = 0; i < 10000; ++i) {
            double log_mid = 0.5 * (log_low + log_high);
            double mid = exp(log_mid);
            calc_tails(a, b, c, d, mid, &lt, &ut);
            double err = fabs(lt - target_alpha);
            if (err < best_err) { best_err = err; best = mid; }
            if (lt > target_alpha) log_low = log_mid; else log_high = log_mid;
            if (log_high - log_low < 1e-32) break;
        }
        *ci_high = best;
    }
}
/* Two-sided exact p-value */
static double exact_p_value(size_t a, size_t b, size_t c, size_t d) {
    size_t r1 = a + b, r2 = c + d, c1 = a + c;
    size_t min_x = (r2 > c1) ? 0 : c1 - r2;
    size_t max_x = (r1 < c1) ? r1 : c1;
    
    double p_obs = exp(log_choose(r1, a) + log_choose(r2, c) - log_choose(r1 + r2, c1));
    double p_val = 0.0;
    
    for (size_t i = min_x; i <= max_x; i++) {
        double p_cur = exp(log_choose(r1, i) + log_choose(r2, c1 - i) - log_choose(r1 + r2, c1));
        if (p_cur <= p_obs * (1.0 + 1e-7)) {
            p_val += p_cur;
        }
    }
    return (p_val > 1.0) ? 1.0 : p_val;
}
/* -----------------------------------------------------------------------
 * Helpers for lm() Linear Regression: OLS Matrix Math & Formula Parsing
 * ----------------------------------------------------------------------- */

/* Gauss-Jordan Elimination with partial pivoting for matrix inversion.
 * Overwrites matrix A with its inverse. Returns 1 on singular failure. */
static int invert_matrix_gj(double *A, UV n) {
	UV *restrict indxc, *restrict indxr, *restrict ipiv;
	Newx(indxc, n, UV);
	Newx(indxr, n, UV);
	Newx(ipiv,  n, UV);

	for (UV i = 0; i < n; i++) ipiv[i] = 0;

	UV icol = 0, irow = 0, ll;
	double big, dum, pivinv, temp;
	bool fail = 0;

	for (UV i = 0; i < n; i++) {
	  big = 0.0;
	  for (UV j = 0; j < n; j++) {
		   if (ipiv[j] != 1) {
		       for (UV k = 0; k < n; k++) {
		           if (ipiv[k] == 0) {
		               if (fabs(A[j * n + k]) >= big) {
		                   big = fabs(A[j * n + k]);
		                   irow = j;
		                   icol = k;
		               }
		           }
		       }
		   }
	  }
	  ++(ipiv[icol]);
	  if (irow != icol) {
		   for (UV l = 0; l < n; l++) {
		       temp = A[irow * n + l];
		       A[irow * n + l] = A[icol * n + l];
		       A[icol * n + l] = temp;
		   }
	  }
	  indxr[i] = irow;
	  indxc[i] = icol;
	  if (A[icol * n + icol] == 0.0) { fail = 1; break; }
	  pivinv = 1.0 / A[icol * n + icol];
	  A[icol * n + icol] = 1.0;
	  for (UV l = 0; l < n; l++) A[icol * n + l] *= pivinv;
	  for (ll = 0; ll < n; ll++) {
		   if (ll != icol) {
		       dum = A[ll * n + icol];
		       A[ll * n + icol] = 0.0;
		       for (UV l = 0; l < n; l++) A[ll * n + l] -= A[icol * n + l] * dum;
		   }
	  }
	}
	if (!fail) {
	  for (UV l = n; l > 0; l--) {
		   if (indxr[l - 1] != indxc[l - 1]) {
		       for (UV k = 0; k < n; k++) {
		           temp = A[k * n + indxr[l - 1]];
		           A[k * n + indxr[l - 1]] = A[k * n + indxc[l - 1]];
		           A[k * n + indxc[l - 1]] = temp;
		       }
		   }
	  }
	}
	Safefree(indxc); Safefree(indxr); Safefree(ipiv);
	return fail;
}

/* Internal extractor resolving single data values from HoA, AoH, or HoH */
static double get_data_value(HV *restrict data_hoa, HV **restrict row_hashes, unsigned int i, const char *restrict var) {
	if (row_hashes) {
	  SV**restrict val = hv_fetch(row_hashes[i], var, strlen(var), 0);
	  if (val && SvROK(*val) && SvTYPE(SvRV(*val)) == SVt_PVAV) {
		   /* Support JSON-decoded HoH where values are arrays: {"mpg":[21.0]} */
		   AV*restrict av = (AV*)SvRV(*val);
		   SV**restrict inner = av_fetch(av, 0, 0);
		   return (inner && SvOK(*inner)) ? SvNV(*inner) : 0.0;
	  }
	  return (val && SvOK(*val)) ? SvNV(*val) : 0.0;
	} else if (data_hoa) {
	  /* Support Hash of Arrays */
	  SV**restrict col = hv_fetch(data_hoa, var, strlen(var), 0);
	  if (col && SvROK(*col) && SvTYPE(SvRV(*col)) == SVt_PVAV) {
		   AV*restrict av = (AV*)SvRV(*col);
		   SV**restrict val = av_fetch(av, i, 0);
		   return (val && SvOK(*val)) ? SvNV(*val) : 0.0;
	  }
	}
	return 0.0;
}

/* Recursive formula resolver evaluating operations for a specific row index */
static double evaluate_term(HV *restrict data_hoa, HV **restrict row_hashes, unsigned int i, const char *restrict term) {
	char term_cpy[256];
	strncpy(term_cpy, term, 255); term_cpy[255] = '\0';

	char *restrict colon = strchr(term_cpy, ':');
	if (colon) {
	  *colon = '\0';
	  return evaluate_term(data_hoa, row_hashes, i, term_cpy) * evaluate_term(data_hoa, row_hashes, i, colon + 1);
	}

	if (strncmp(term_cpy, "I(", 2) == 0) {
	  char *restrict end = strrchr(term_cpy, ')');
	  if (end) *end = '\0';
	  char *restrict inner = term_cpy + 2;
	  char *restrict caret = strchr(inner, '^');
	  int power = 1;
	  if (caret) {
		   *caret = '\0';
		   power = atoi(caret + 1);
	  }
	  double v = get_data_value(data_hoa, row_hashes, i, inner);
	  return power == 1 ? v : pow(v, power);
	}
	return get_data_value(data_hoa, row_hashes, i, term_cpy);
}

/* Struct for sorting p-values while remembering their original index */
typedef struct {
	double p;
	size_t orig_idx;
} PVal;

/* Comparator for qsort */
static int cmp_pval(const void *a, const void *b) {
    double diff = ((PVal*)a)->p - ((PVal*)b)->p;
    if (diff < 0) return -1;
    if (diff > 0) return 1;
    /* Stabilize sort by falling back to original index */
    return ((PVal*)a)->orig_idx - ((PVal*)b)->orig_idx; 
}
/* -----------------------------------------------------------------------
 * Helpers for cor(): ranking (Spearman), Pearson r, Kendall tau-b
 * ----------------------------------------------------------------------- */
/* Item used to sort values while remembering their original index,
 * needed for average-rank tie-breaking in Spearman correlation.        */
typedef struct {
    double val;
    size_t idx;
} RankItem;

static int cmp_rank_item(const void *a, const void *b) {
    double diff = ((RankItem*)a)->val - ((RankItem*)b)->val;
    if (diff < 0) return -1;
    if (diff > 0) return  1;
    return 0;
}

/* Compute 1-based average ranks with tie-breaking into out[].
 * in[] is not modified.                                                 */
static void rank_data(const double *in, double *out, size_t n) {
	RankItem *restrict ri;
	Newx(ri, n, RankItem);
	for (size_t i = 0; i < n; i++) { ri[i].val = in[i]; ri[i].idx = i; }
	qsort(ri, n, sizeof(RankItem), cmp_rank_item);

	size_t i = 0;
	while (i < n) {
	  size_t j = i;
	  /* Find the full extent of this tie group */
	  while (j + 1 < n && ri[j + 1].val == ri[j].val) j++;
	  /* All members get the average of ranks i+1 … j+1 (1-based) */
	  double avg = (double)(i + j) / 2.0 + 1.0;
	  for (size_t k = i; k <= j; k++) out[ri[k].idx] = avg;
	  i = j + 1;
	}
	Safefree(ri);
}

/* Pearson product-moment r between two n-element arrays.
 * Returns NAN when either variable has zero variance (matches R).       */
static double pearson_corr(const double *restrict x, const double *restrict y, unsigned int n) {
    double sx = 0, sy = 0, sxy = 0, sx2 = 0, sy2 = 0;
    for (size_t i = 0; i < n; i++) {
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
    size_t C = 0, D = 0, tie_x = 0, tie_y = 0;
    for (size_t i = 0; i < n - 1; i++) {
        for (size_t j = i + 1; j < n; j++) {
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
                           size_t n, const char *method) {
	if (strcmp(method, "spearman") == 0) {
	  double *restrict rx, *restrict ry;
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

int compare_doubles(const void *a, const void *b) {
    double da = *(const double*restrict)a;
    double db = *(const double*restrict)b;
    return (da > db) - (da < db);
}
/* --- XS SECTION --- */
MODULE = Stats::LikeR	PACKAGE = Stats::LikeR		

double mean(...)
	PROTOTYPE: @
	INIT:
	  double total = 0;
	  size_t count = 0;
	CODE:
	  for (unsigned int i = 0; i < items; i++) {
		   SV*restrict arg = ST(i);
		   if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		       AV*restrict av = (AV*)SvRV(arg);
		       size_t len = av_len(av) + 1;
		       for (size_t j = 0; j < len; j++) {
		           SV**restrict tv = av_fetch(av, j, 0);
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

double sd(...)
	PROTOTYPE: @
	INIT:
		double total = 0, sum_sq = 0;
		size_t count = 0;
	CODE:
	  for (size_t i = 0; i < items; i++) {
		   SV* arg = ST(i);
		   if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		       AV* av = (AV*)SvRV(arg);
		       size_t len = av_len(av) + 1;
		       for (size_t j = 0; j < len; j++) {
		           SV**restrict tv = av_fetch(av, j, 0);
		           if (tv && SvOK(*tv)) { total += SvNV(*tv); count++; }
		       }
		   } else if (SvOK(arg)) {
		       total += SvNV(arg); count++;
		   }
	  }
	  if (count < 2) croak("stdev needs >= 2 elements");
	  double avg = total / count;
	  for (size_t i = 0; i < items; i++) {
		   SV*restrict arg = ST(i);
		   if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		       AV*restrict av = (AV*)SvRV(arg);
		       size_t len = av_len(av) + 1;
		       for (size_t j = 0; j < len; j++) {
		           SV**restrict tv = av_fetch(av, j, 0);
		           if (tv && SvOK(*tv)) sum_sq += pow(SvNV(*tv) - avg, 2);
		       }
		   } else if (SvOK(arg)) {
		       sum_sq += pow(SvNV(arg) - avg, 2);
		   }
	  }
	  RETVAL = sqrt(sum_sq / (count - 1));
	OUTPUT:
	  RETVAL

double var(...)
	PROTOTYPE: @
	INIT:
	  double total = 0.0, sum_sq = 0.0;
	  size_t count = 0;
	CODE:
	  for (size_t i = 0; i < items; i++) {
		   SV*restrict arg = ST(i);
		   if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		       AV*restrict av = (AV*)SvRV(arg);
		       size_t len = av_len(av) + 1;
		       for (size_t j = 0; j < len; j++) {
		           SV**restrict tv = av_fetch(av, j, 0);
		           if (tv && SvOK(*tv)) { total += SvNV(*tv); count++; }
		       }
		   } else if (SvOK(arg)) {
		       total += SvNV(arg); count++;
		   }
	  }
	  if (count < 2) croak("stdev needs >= 2 elements");
	  double avg = total / count;
	  for (size_t i = 0; i < items; i++) {
		   SV*restrict arg = ST(i);
		   if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		       AV*restrict av = (AV*)SvRV(arg);
		       size_t len = av_len(av) + 1;
		       for (ssize_t j = 0; j < len; j++) {
		           SV**restrict tv = av_fetch(av, j, 0);
		           if (tv && SvOK(*tv)) sum_sq += pow(SvNV(*tv) - avg, 2);
		       }
		   } else if (SvOK(arg)) {
		       sum_sq += pow(SvNV(arg) - avg, 2);
		   }
	  }
	  RETVAL = sum_sq / (count - 1);
	OUTPUT:
	  RETVAL

double pearson_r(SV* x_sv, SV* y_sv)
	INIT:
	  if (!SvROK(x_sv) || SvTYPE(SvRV(x_sv)) != SVt_PVAV ||
		   !SvROK(y_sv) || SvTYPE(SvRV(y_sv)) != SVt_PVAV) {
		   croak("Both arguments must be array references");
	  }
	  AV*restrict x_av = (AV*)SvRV(x_sv);
	  AV*restrict y_av = (AV*)SvRV(y_sv);
	  size_t n = av_len(x_av) + 1;
	  if (n != (av_len(y_av) + 1)) croak("Arrays must have the same number of elements");
	  if (n < 2) croak("Need at least 2 elements");
	  double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0, sum_y2 = 0;
	CODE:
	  for (size_t i = 0; i < n; i++) {
		   SV**restrict xv = av_fetch(x_av, i, 0);
		   SV**restrict yv = av_fetch(y_av, i, 0);
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


SV* t_test(...)
	CODE:
	{
	  if (items % 2 != 0)
		   croak("Usage: t_test(x => [...], y => [...], ...) - must be even key/value pairs");

	  /* --- Parse named arguments from the flat stack --- */
	  SV*restrict x_sv = NULL;
	  SV*restrict y_sv = NULL;
	  double mu = 0.0, conf_level = 0.95;
	  bool paired = FALSE, var_equal = FALSE;
	  const char*restrict alternative = "two.sided";

	  for (I32 i = 0; i < items; i += 2) {
		   const char*restrict key = SvPV_nolen(ST(i));
		   SV*restrict val = ST(i + 1);

		   if      (strEQ(key, "x"))           x_sv        = val;
		   else if (strEQ(key, "y"))           y_sv        = val;
		   else if (strEQ(key, "mu"))          mu          = SvNV(val);
		   else if (strEQ(key, "paired"))      paired      = SvTRUE(val);
		   else if (strEQ(key, "var_equal"))   var_equal   = SvTRUE(val);
		   else if (strEQ(key, "conf_level"))  conf_level  = SvNV(val);
		   else if (strEQ(key, "alternative")) alternative = SvPV_nolen(val);
		   else croak("t_test: unknown argument '%s'", key);
	  }

	  /* --- Validate required / types --- */
	  if (!x_sv || !SvROK(x_sv) || SvTYPE(SvRV(x_sv)) != SVt_PVAV)
		   croak("t_test: 'x' is a required argument and must be an ARRAY reference");

	  AV*restrict x_av = (AV*)SvRV(x_sv);
	  size_t nx = av_len(x_av) + 1;
	  if (nx < 2) croak("t_test: 'x' needs at least 2 elements");

	  AV*restrict y_av = NULL;
	  if (y_sv && SvROK(y_sv) && SvTYPE(SvRV(y_sv)) == SVt_PVAV)
		   y_av = (AV*)SvRV(y_sv);

	  if (conf_level <= 0.0 || conf_level >= 1.0)
		   croak("t_test: 'conf_level' must be between 0 and 1");

	  /* --- Computation (identical to the original) --- */
	  double sum_x = 0, sum_x2 = 0, mean_x, var_x;
	  double t_stat, df, p_val, std_err, cint_est;
	  HV*restrict results = newHV();

	  for (size_t i = 0; i < nx; i++) {
		   SV**restrict tv = av_fetch(x_av, i, 0);
		   double val = (tv && SvOK(*tv)) ? SvNV(*tv) : 0;
		   sum_x += val; sum_x2 += val * val;
	  }
	  mean_x = sum_x / nx;
	  var_x  = (sum_x2 - (sum_x * sum_x) / nx) / (nx - 1);

	  if (paired || y_av) {
		   if (!y_av) croak("t_test: 'y' must be provided for paired or two-sample tests");
		   size_t ny = av_len(y_av) + 1;
		   if (paired && ny != nx) croak("t_test: Paired arrays must be same length");

		   double sum_y = 0, sum_y2 = 0, mean_y, var_y;
		   for (size_t i = 0; i < ny; i++) {
		       SV**restrict tv = av_fetch(y_av, i, 0);
		       double val = (tv && SvOK(*tv)) ? SvNV(*tv) : 0;
		       sum_y += val; sum_y2 += val * val;
		   }
		   mean_y = sum_y / ny;
		   var_y  = (sum_y2 - (sum_y * sum_y) / ny) / (ny - 1);

		   if (paired) {
		       double sum_d = 0, sum_d2 = 0;
		       for (size_t i = 0; i < nx; i++) {
		           double dx = SvNV(*av_fetch(x_av, i, 0));
		           double dy = SvNV(*av_fetch(y_av, i, 0));
		           sum_d += (dx - dy); sum_d2 += (dx - dy) * (dx - dy);
		       }
		       double mean_d = sum_d / nx;
		       double var_d  = (sum_d2 - (sum_d * sum_d) / nx) / (nx - 1);
		       cint_est = mean_d;
		       std_err  = sqrt(var_d / nx);
		       t_stat   = (cint_est - mu) / std_err;
		       df       = nx - 1;
		       hv_store(results, "estimate", 8, newSVnv(mean_d), 0);
		   } else if (var_equal) {
		       double pooled_var = ((nx - 1) * var_x + (ny - 1) * var_y) / (nx + ny - 2);
		       cint_est = mean_x - mean_y;
		       std_err  = sqrt(pooled_var * (1.0 / nx + 1.0 / ny));
		       t_stat   = (cint_est - mu) / std_err;
		       df       = nx + ny - 2;
		       hv_store(results, "estimate_x", 10, newSVnv(mean_x), 0);
		       hv_store(results, "estimate_y", 10, newSVnv(mean_y), 0);
		   } else {
		       cint_est         = mean_x - mean_y;
		       double stderr_x2 = var_x / nx;
		       double stderr_y2 = var_y / ny;
		       std_err          = sqrt(stderr_x2 + stderr_y2);
		       t_stat           = (cint_est - mu) / std_err;
		       df = pow(stderr_x2 + stderr_y2, 2) /
		            (pow(stderr_x2, 2) / (nx - 1) + pow(stderr_y2, 2) / (ny - 1));
		       hv_store(results, "estimate_x", 10, newSVnv(mean_x), 0);
		       hv_store(results, "estimate_y", 10, newSVnv(mean_y), 0);
		   }
	  } else {
		   cint_est = mean_x;
		   std_err  = sqrt(var_x / nx);
		   t_stat   = (cint_est - mu) / std_err;
		   df       = nx - 1;
		   hv_store(results, "estimate", 8, newSVnv(mean_x), 0);
	  }

	  p_val = get_t_pvalue(t_stat, df, alternative);

	  double alpha = 1.0 - conf_level;
	  double t_crit, ci_lower, ci_upper;

	  if (strcmp(alternative, "less") == 0) {
		   t_crit   = qt_tail(df, alpha);
		   ci_lower = -INFINITY;
		   ci_upper = cint_est + t_crit * std_err;
	  } else if (strcmp(alternative, "greater") == 0) {
		   t_crit   = qt_tail(df, alpha);
		   ci_lower = cint_est - t_crit * std_err;
		   ci_upper = INFINITY;
	  } else {
		   t_crit   = qt_tail(df, alpha / 2.0);
		   ci_lower = cint_est - t_crit * std_err;
		   ci_upper = cint_est + t_crit * std_err;
	  }

	  AV*restrict conf_int = newAV();
	  av_push(conf_int, newSVnv(ci_lower));
	  av_push(conf_int, newSVnv(ci_upper));

	  hv_store(results, "statistic", 9, newSVnv(t_stat), 0);
	  hv_store(results, "df",        2, newSVnv(df),     0);
	  hv_store(results, "p_value",   7, newSVnv(p_val),  0);
	  hv_store(results, "conf_int",  8, newRV_noinc((SV*)conf_int), 0);

	  RETVAL = newRV_noinc((SV*)results);
	}
	OUTPUT:
	  RETVAL

void p_adjust(SV* p_sv, const char* method = "holm")
	INIT:
	  if (!SvROK(p_sv) || SvTYPE(SvRV(p_sv)) != SVt_PVAV) {
		   croak("p_adjust: first argument must be an ARRAY reference of p-values");
	  }
	  AV *restrict p_av = (AV*)SvRV(p_sv);
	  size_t n = av_len(p_av) + 1;
	  /* Handle empty input */
	  if (n == 0) {
		   XSRETURN_EMPTY;
	  }
	  /* Normalize method string */
	  char meth[64];
	  strncpy(meth, method, 63); meth[63] = '\0';
	  for(unsigned short int i = 0; meth[i]; i++) meth[i] = tolower(meth[i]);
	  /* Resolve aliases */
	  if (strstr(meth, "benjamini") && strstr(meth, "hochberg")) strcpy(meth, "bh");
	  if (strstr(meth, "benjamini") && strstr(meth, "yekutieli")) strcpy(meth, "by");
	  if (strcmp(meth, "fdr") == 0) strcpy(meth, "bh");
	  /* Allocate C memory */
	  PVal *restrict arr;
	  double *restrict adj;
	  Newx(arr, n, PVal);
	  Newx(adj, n, double);

	  for (size_t i = 0; i < n; i++) {
		   SV**restrict tv = av_fetch(p_av, i, 0);
		   arr[i].p = (tv && SvOK(*tv)) ? SvNV(*tv) : 1.0;
		   arr[i].orig_idx = i;
	  }

	  /* Sort ascending (Stable sort using original index) */
	  qsort(arr, n, sizeof(PVal), cmp_pval);

	PPCODE:
	  if (strcmp(meth, "bonferroni") == 0) {
		   for (size_t i = 0; i < n; i++) {
		       double v = arr[i].p * n;
		       adj[arr[i].orig_idx] = (v < 1.0) ? v : 1.0;
		   }
	  } else if (strcmp(meth, "holm") == 0) {
		   double cummax = 0.0;
		   for (UV i = 0; i < n; i++) {
		       double v = arr[i].p * (n - i);
		       if (v > cummax) cummax = v;
		       adj[arr[i].orig_idx] = (cummax < 1.0) ? cummax : 1.0;
		   }
	  } else if (strcmp(meth, "hochberg") == 0) {
		   double cummin = 1.0;
		   for (ssize_t i = n - 1; i >= 0; i--) {
		       double v = arr[i].p * (n - i);
		       if (v < cummin) cummin = v;
		       adj[arr[i].orig_idx] = (cummin < 1.0) ? cummin : 1.0;
		   }
	  } else if (strcmp(meth, "bh") == 0) {
		   double cummin = 1.0;
		   for (ssize_t i = n - 1; i >= 0; i--) {
		       double v = arr[i].p * n / (i + 1.0);
		       if (v < cummin) cummin = v;
		       adj[arr[i].orig_idx] = (cummin < 1.0) ? cummin : 1.0;
		   }
	  } else if (strcmp(meth, "by") == 0) {
		   double q = 0.0;
		   for (size_t i = 1; i <= n; i++) q += 1.0 / i;
		   double cummin = 1.0;
		   for (ssize_t i = n - 1; i >= 0; i--) {
		       double v = arr[i].p * n / (i + 1.0) * q;
		       if (v < cummin) cummin = v;
		       adj[arr[i].orig_idx] = (cummin < 1.0) ? cummin : 1.0;
		   }
	  } else if (strcmp(meth, "hommel") == 0) {
		   double *restrict pa = (double *)malloc(n * sizeof(double));
		   double *restrict q_arr = (double *)malloc(n * sizeof(double));
		   
		   /* Initial: min(n * p[i] / (i + 1)) */
		   double min_val = n * arr[0].p;
		   for (size_t i = 1; i < n; i++) {
		       double temp = (n * arr[i].p) / (i + 1.0);
		       if (temp < min_val) {
		           min_val = temp;
		       }
		   }
		   /* pa <- q <- rep(min, n) */
		   for (size_t i = 0; i < n; i++) {
		       pa[i] = min_val;
		       q_arr[i] = min_val;
		   }
		   for (size_t j = n - 1; j >= 2; j--) {
		       ssize_t n_mj = n - j;       /* Max index for 'ij'. Length is n_mj + 1 */
		       ssize_t i2_len = j - 1;     /* Length of 'i2' */
		       /* Calculate q1 = min(j * p[i2] / (2:j)) */
		       double q1 = (j * arr[n_mj + 1].p) / 2.0;
		       for (size_t k = 1; k < i2_len; k++) {
		           double temp_q1 = (j * arr[n_mj + 1 + k].p) / (2.0 + k);
		           if (temp_q1 < q1) {
		               q1 = temp_q1;
		           }
		       }
		       /* q[ij] <- pmin(j * p[ij], q1) */
		       for (size_t i = 0; i <= n_mj; i++) {
		           double v = j * arr[i].p;
		           q_arr[i] = (v < q1) ? v : q1;
		       }
		       /* q[i2] <- q[n - j] */
		       for (size_t i = 0; i < i2_len; i++) {
		           q_arr[n_mj + 1 + i] = q_arr[n_mj];
		       }
		       /* pa <- pmax(pa, q) */
		       for (size_t i = 0; i < n; i++) {
		           if (pa[i] < q_arr[i]) {
		              pa[i] = q_arr[i];
		           }
		       }
		   }
		   
		   /* pmin(1, pmax(pa, p))[ro] — map sorted results back to original indices */
		   for (size_t i = 0; i < n; i++) {
		       double v = (pa[i] > arr[i].p) ? pa[i] : arr[i].p;
		       if (v > 1.0) v = 1.0;
		       adj[arr[i].orig_idx] = v;
		   }
		   free(pa); pa = NULL;
		   free(q_arr); q_arr = NULL;
	  } else if (strcmp(meth, "none") == 0) {
		   for (size_t i = 0; i < n; i++) {
		   	adj[arr[i].orig_idx] = arr[i].p;
		   }
	  } else {
		   Safefree(arr); Safefree(adj);
		   croak("Unknown p-value adjustment method: %s", method);
	  }
	  /* Push values onto the Perl stack as a flat list */
	  EXTEND(SP, n);
	  for (size_t i = 0; i < n; i++) {
		   PUSHs(sv_2mortal(newSVnv(adj[i])));
	  }
	  Safefree(arr); arr = NULL;
	  Safefree(adj); adj = NULL;

double median(...)
	PROTOTYPE: @
	INIT:
	  size_t total_count = 0, k = 0;
	  double *restrict nums;
	  double median_val = 0.0;
	CODE:
	  /* Pass 1: Count total valid elements to allocate memory */
	  for (size_t i = 0; i < items; i++) {
		   SV*restrict arg = ST(i);
		   if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		       AV* av = (AV*)SvRV(arg);
		       unsigned int len = av_len(av) + 1;
		       for (size_t j = 0; j < len; j++) {
		           SV**restrict tv = av_fetch(av, j, 0);
		           if (tv && SvOK(*tv)) { total_count++; }
		       }
		   } else if (SvOK(arg)) {
		       total_count++;
		   }
	  }
	  if (total_count == 0) croak("median needs >= 1 element");
	  /* Allocate C array now that we know the exact size */
	  Newx(nums, total_count, double);
	  /* Pass 2: Populate the C array */
	  for (size_t i = 0; i < items; i++) {
		   SV*restrict arg = ST(i);
		   if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		       AV*restrict av = (AV*)SvRV(arg);
		       size_t len = av_len(av) + 1;
		       for (size_t j = 0; j < len; j++) {
		           SV**restrict tv = av_fetch(av, j, 0);
		           if (tv && SvOK(*tv)) { 
		               nums[k++] = SvNV(*tv); 
		           }
		       }
		   } else if (SvOK(arg)) {
		       nums[k++] = SvNV(arg);
		   }
	  }
	  /* Sort and calculate median */
	  qsort(nums, total_count, sizeof(double), compare_doubles);
	  if (total_count % 2 == 0) {
		   median_val = (nums[total_count / 2 - 1] + nums[total_count / 2]) / 2.0;
	  } else {
		   median_val = nums[total_count / 2];
	  }
	  Safefree(nums); nums = NULL;
	  RETVAL = median_val;
	OUTPUT:
	  RETVAL

SV* cor(SV* x_sv, SV* y_sv = &PL_sv_undef, const char* method = "pearson")
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
	  AV*restrict x_av = (AV*)SvRV(x_sv);
	  size_t nx   = av_len(x_av) + 1;
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
	  AV*restrict y_av     = has_y ? (AV*)SvRV(y_sv) : NULL;
	  size_t ny       = has_y ? av_len(y_av) + 1 : 0;
	  short int y_is_matrix = 0;
	  if (has_y && ny > 0) {
		   SV** fp = av_fetch(y_av, 0, 0);
		   if (fp && SvROK(*fp) && SvTYPE(SvRV(*fp)) == SVt_PVAV)
		       y_is_matrix = 1;
	  }
	CODE:// Branch 1: both inputs are flat vectors  →  scalar result
	  if (!x_is_matrix && !y_is_matrix) {
		   if (!has_y) {
		       /* cor(vector) == 1 by definition */
		       RETVAL = newSVnv(1.0);
		   } else {
		       if (nx != ny)
		           croak("cor: x and y must have the same length (%lu vs %lu)",
		                 nx, ny);
		       if (nx < 2)
		           croak("cor: need at least 2 observations");

		       double *restrict xd, *restrict yd;
		       Newx(xd, nx, double);
		       Newx(yd, ny, double);

		       for (size_t i = 0; i < nx; i++) {
		           SV**restrict tv = av_fetch(x_av, i, 0);
		           xd[i] = (tv && SvOK(*tv)) ? SvNV(*tv) : 0.0;
		       }
		       for (size_t i = 0; i < ny; i++) {
		           SV**restrict tv = av_fetch(y_av, i, 0);
		           yd[i] = (tv && SvOK(*tv)) ? SvNV(*tv) : 0.0;
		       }

		       double r = compute_cor(xd, yd, nx, method);
		       Safefree(xd); Safefree(yd);
		       RETVAL = newSVnv(r);
		   }
	  } else {//Branch 2: x is a matrix (or y is a matrix)  →  AoA result
		   /* -- resolve x matrix dimensions ---------------------------- */
		   if (!x_is_matrix)
		       croak("cor: x must be a matrix (array ref of array refs) "
		             "when y is a matrix");

		   SV**restrict xr0 = av_fetch(x_av, 0, 0);
		   if (!xr0 || !SvROK(*xr0) || SvTYPE(SvRV(*xr0)) != SVt_PVAV)
		       croak("cor: each row of x must be an ARRAY reference");
		   size_t ncols_x = av_len((AV*)SvRV(*xr0)) + 1;
		   if (ncols_x == 0) croak("cor: x matrix has zero columns");
		   size_t nrows   = nx;   /* observations */

		   /* -- extract x columns -------------------------------------- */
		   double **restrict col_x;
		   Newx(col_x, ncols_x, double*);
		   for (size_t j = 0; j < ncols_x; j++) {
		       Newx(col_x[j], nrows, double);
		       for (size_t i = 0; i < nrows; i++) {
		           SV**restrict rv = av_fetch(x_av, i, 0);
		           if (!rv || !SvROK(*rv))
		               croak("cor: x row %lu is not an array ref", i);
		           AV*restrict  row = (AV*)SvRV(*rv);
		           SV**restrict cv  = av_fetch(row, j, 0);
		           col_x[j][i] = (cv && SvOK(*cv)) ? SvNV(*cv) : 0.0;
		       }
		   }

		   /* -- resolve y: separate matrix or re-use x (symmetric) ---- */
		   size_t    ncols_y;
		   double **restrict col_y   = NULL;
		   int    symmetric = 0;   /* 1 = cor(X) — result is symmetric */
		   if (has_y && y_is_matrix) {
		       /* cross-correlation: X (nrows × p) vs Y (nrows × q)      */
		       if (ny != nrows)
		           croak("cor: x and y must have the same number of rows "
		                 "(%lu vs %lu)", nrows, ny);
		       SV**restrict yr0 = av_fetch(y_av, 0, 0);
		       if (!yr0 || !SvROK(*yr0) || SvTYPE(SvRV(*yr0)) != SVt_PVAV)
		           croak("cor: each row of y must be an ARRAY reference");
		       ncols_y = av_len((AV*)SvRV(*yr0)) + 1;
		       if (ncols_y == 0) croak("cor: y matrix has zero columns");

		       Newx(col_y, ncols_y, double*);
		       for (unsigned int j = 0; j < ncols_y; j++) {
		           Newx(col_y[j], nrows, double);
		           for (unsigned int i = 0; i < nrows; i++) {
		               SV**  rv = av_fetch(y_av, i, 0);
		               if (!rv || !SvROK(*rv))
		                   croak("cor: y row %d is not an array ref", i);
		               AV*restrict  row = (AV*)SvRV(*rv);
		               SV**restrict cv  = av_fetch(row, j, 0);
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
		       croak("cor: need at least 2 observations (got %lu)", nrows);
		   /* -- build cache for symmetric case: compute upper triangle,
		      store results, mirror to lower triangle                    */
		   AV*restrict result_av = newAV();
		   av_extend(result_av, ncols_x - 1);
		   /* Allocate per-row AVs up front so we can fill them in order */
		   AV **restrict rows_out;
		   Newx(rows_out, ncols_x, AV*);
		   for (size_t i = 0; i < ncols_x; i++) {
		       rows_out[i] = newAV();
		       av_extend(rows_out[i], ncols_y - 1);
		   }
		   if (symmetric) {
		       /* Upper triangle + diagonal, then mirror.
		          r_cache[i][j] (j >= i) holds the computed value.      */
		       double **restrict r_cache;
		       Newx(r_cache, ncols_x, double*);
		       for (size_t i = 0; i < ncols_x; i++)
		           Newx(r_cache[i], ncols_x, double);

		       for (size_t i = 0; i < ncols_x; i++) {
		           r_cache[i][i] = 1.0;                /* diagonal */
		           for (size_t j = i + 1; j < ncols_x; j++) {
		               double r = compute_cor(col_x[i], col_x[j], nrows, method);
		               r_cache[i][j] = r;
		               r_cache[j][i] = r;              /* symmetry */
		           }
		       }
		       /* fill output AoA from cache */
		       for (size_t i = 0; i < ncols_x; i++)
		           for (size_t j = 0; j < ncols_x; j++)
		               av_store(rows_out[i], j, newSVnv(r_cache[i][j]));

		       for (size_t i = 0; i < ncols_x; i++) Safefree(r_cache[i]);
		       Safefree(r_cache);
		   } else {
		       /* cross-correlation: every (i,j) pair is independent     */
		       for (size_t i = 0; i < ncols_x; i++)
		           for (size_t j = 0; j < ncols_y; j++)
		               av_store(rows_out[i], j,
		                        newSVnv(compute_cor(col_x[i], col_y[j],
		                                            nrows, method)));
		   }
		   /* push row AVs into result */
		   for (size_t i = 0; i < ncols_x; i++)
		       av_store(result_av, i, newRV_noinc((SV*)rows_out[i]));
		   Safefree(rows_out);
		   /* -- free column arrays ------------------------------------- */
		   for (size_t j = 0; j < ncols_x; j++) Safefree(col_x[j]);
		   Safefree(col_x);
		   if (!symmetric) {
		       for (size_t j = 0; j < ncols_y; j++) Safefree(col_y[j]);
		       Safefree(col_y);
		   }
		   RETVAL = newRV_noinc((SV*)result_av);
	  }
	OUTPUT:
	  RETVAL

void scale(...)
	PROTOTYPE: @
	PPCODE:
	{
	  bool do_center_mean = true, do_scale_sd = true;
	  double center_val = 0.0, scale_val = 1.0;
	  int data_items = items;
	  /* 1. Parse Options Hash (if it exists as the last argument) */
	  if (items > 0) {
		   SV*restrict last_arg = ST(items - 1);
		   if (SvROK(last_arg) && SvTYPE(SvRV(last_arg)) == SVt_PVHV) {
		       data_items = items - 1; /* Exclude hash from data processing */
		       HV* opt_hv = (HV*)SvRV(last_arg);

		       /* --- Parse 'center' --- */
		       SV**restrict center_sv = hv_fetch(opt_hv, "center", 6, 0);
		       if (center_sv) {
		           SV*restrict val_sv = *center_sv;
		           if (!SvOK(val_sv)) {
		               do_center_mean = false; center_val = 0.0;
		           } else {
		               char *restrict str = SvPV_nolen(val_sv);
		               /* Trap booleans and empty strings before numeric checks */
		               if (strcasecmp(str, "mean") == 0 || strcasecmp(str, "true") == 0 || strcmp(str, "1") == 0) {
		                   do_center_mean = true;
		               } else if (strcasecmp(str, "none") == 0 || strcasecmp(str, "false") == 0 || strcmp(str, "0") == 0 || strcmp(str, "") == 0) {
		                   do_center_mean = false; center_val = 0.0;
		               } else if (looks_like_number(val_sv)) {
		                   do_center_mean = false; center_val = SvNV(val_sv);
		               } else if (SvTRUE(val_sv)) {
		                   do_center_mean = true;
		               } else {
		                   do_center_mean = false; center_val = 0.0;
		               }
		           }
		       }
		       /* --- Parse 'scale' --- */
		       SV**restrict scale_sv = hv_fetch(opt_hv, "scale", 5, 0);
		       if (scale_sv) {
		           SV*restrict val_sv = *scale_sv;
		           if (!SvOK(val_sv)) {
		               do_scale_sd = false; scale_val = 1.0;
		           } else {
		               char *str = SvPV_nolen(val_sv);
		               if (strcasecmp(str, "sd") == 0 || strcasecmp(str, "true") == 0 || strcmp(str, "1") == 0) {
		                   do_scale_sd = true;
		               } else if (strcasecmp(str, "none") == 0 || strcasecmp(str, "false") == 0 || strcmp(str, "0") == 0 || strcmp(str, "") == 0) {
		                   do_scale_sd = false; scale_val = 1.0;
		               } else if (looks_like_number(val_sv)) {
		                   do_scale_sd = false; scale_val = SvNV(val_sv);
		                   if (scale_val == 0.0) scale_val = 1.0; /* Prevent Division By Zero */
		               } else if (SvTRUE(val_sv)) {
		                   do_scale_sd = true;
		               } else {
		                   do_scale_sd = false; scale_val = 1.0;
		               }
		           }
		       }
		   }
	  }
	  /* 2. Detect if the input is a Matrix (Array of Arrays) */
	  bool is_matrix = false;
	  if (data_items == 1) {
		   SV*restrict first_arg = ST(0);
		   if (SvROK(first_arg) && SvTYPE(SvRV(first_arg)) == SVt_PVAV) {
		       AV*restrict av = (AV*)SvRV(first_arg);
		       if (av_len(av) >= 0) {
		           SV**restrict first_elem = av_fetch(av, 0, 0);
		           if (first_elem && SvROK(*first_elem) && SvTYPE(SvRV(*first_elem)) == SVt_PVAV) {
		               is_matrix = true;
		           }
		       }
		   }
	  }
	  if (is_matrix) {
		   /* ========================================================= */
		   /* MATRIX MODE: Scale columns independently (Just like R)    */
		   /* ========================================================= */
		   AV*restrict mat_av = (AV*)SvRV(ST(0));
		   size_t nrow = av_len(mat_av) + 1, ncol = 0;
		   
		   SV**restrict first_row = av_fetch(mat_av, 0, 0);
		   ncol = av_len((AV*)SvRV(*first_row)) + 1;

		   if (nrow == 0 || ncol == 0) croak("scale requires non-empty matrix");

		   /* Create a new matrix for the scaled output */
		   AV*restrict result_av = newAV();
		   av_extend(result_av, nrow - 1);
		   AV**restrict row_ptrs = (AV**)safemalloc(nrow * sizeof(AV*));
		   for (size_t r = 0; r < nrow; r++) {
		       row_ptrs[r] = newAV();
		       av_extend(row_ptrs[r], ncol - 1);
		       av_push(result_av, newRV_noinc((SV*)row_ptrs[r]));
		   }

		   /* Calculate and apply scale per column */
		   for (size_t c = 0; c < ncol; c++) {
		       double col_sum = 0.0;
		       double *restrict col_data;
		       Newx(col_data, nrow, double);

		       /* Extract the column data */
		       for (size_t r = 0; r < nrow; r++) {
		           SV**restrict row_sv = av_fetch(mat_av, r, 0);
		           if (row_sv && SvROK(*row_sv)) {
		               AV*restrict row_av = (AV*)SvRV(*row_sv);
		               SV**restrict cell_sv = av_fetch(row_av, c, 0);
		               col_data[r] = (cell_sv && SvOK(*cell_sv)) ? SvNV(*cell_sv) : 0.0;
		           } else {
		               col_data[r] = 0.0;
		           }
		           col_sum += col_data[r];
		       }

		       double col_center = do_center_mean ? (col_sum / nrow) : center_val;
		       double col_scale = scale_val;
		       /* Calculate Standard Deviation for this specific column if needed */
		       if (do_scale_sd) {
		           if (nrow <= 1) {
		               Safefree(col_data);
		               safefree(row_ptrs);
		               croak("scale needs >= 2 rows to calculate standard deviation for a matrix column");
		           }
		           double sum_sq = 0.0;
		           for (size_t r = 0; r < nrow; r++) {
		               double diff = col_data[r] - col_center;
		               sum_sq += diff * diff;
		           }
		           col_scale = sqrt(sum_sq / (nrow - 1));
		       }
		       /* Store scaled values back into the new matrix rows */
		       for (size_t r = 0; r < nrow; r++) {
		           double centered = col_data[r] - col_center;
		           double final_val = (col_scale == 0.0) ? (0.0 / 0.0) : (centered / col_scale);
		           av_store(row_ptrs[r], c, newSVnv(final_val));
		       }
		       Safefree(col_data);
		   }
		   safefree(row_ptrs);
		   /* Push the resulting matrix as a single Reference onto the Perl stack */
		   EXTEND(SP, 1);
		   PUSHs(sv_2mortal(newRV_noinc((SV*)result_av)));
	  } else {
		   /* ========================================================= */
		   /* FLAT LIST MODE: Original functionality                    */
		   /* ========================================================= */
		   size_t total_count = 0, k = 0;
		   double *restrict nums;
		   double sum = 0.0;
		   for (size_t i = 0; i < data_items; i++) {
		       SV* arg = ST(i);
		       if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		           AV*restrict av = (AV*)SvRV(arg);
		           size_t len = av_len(av) + 1;
		           for (unsigned int j = 0; j < len; j++) {
		               SV** tv = av_fetch(av, j, 0);
		               if (tv && SvOK(*tv)) { total_count++; }
		           }
		       } else if (SvOK(arg)) {
		           total_count++;
		       }
		   }

		   if (total_count == 0) croak("scale requires at least 1 numeric element");
		   Newx(nums, total_count, double);

		   for (size_t i = 0; i < data_items; i++) {
		       SV*restrict arg = ST(i);
		       if (SvROK(arg) && SvTYPE(SvRV(arg)) == SVt_PVAV) {
		           AV*restrict av = (AV*)SvRV(arg);
		           size_t len = av_len(av) + 1;
		           for (size_t j = 0; j < len; j++) {
		               SV**restrict tv = av_fetch(av, j, 0);
		               if (tv && SvOK(*tv)) { 
		                   double val = SvNV(*tv);
		                   nums[k++] = val; sum += val;
		               }
		           }
		       } else if (SvOK(arg)) {
		           double val = SvNV(arg);
		           nums[k++] = val; sum += val;
		       }
		   }

		   if (do_center_mean) center_val = sum / total_count;

		   if (do_scale_sd) {
		       if (total_count <= 1) {
		           Safefree(nums);
		           croak("scale needs >= 2 elements to calculate SD");
		       }
		       double sum_sq = 0.0;
		       for (size_t i = 0; i < total_count; i++) {
		           double diff = nums[i] - center_val;
		           sum_sq += diff * diff;
		       }
		       scale_val = sqrt(sum_sq / (total_count - 1));
		   }

		   EXTEND(SP, total_count);
		   for (size_t i = 0; i < total_count; i++) {
		       double centered = nums[i] - center_val;
		       double final_val = (scale_val == 0.0) ? (0.0 / 0.0) : (centered / scale_val);
		       PUSHs(sv_2mortal(newSVnv(final_val)));
		   }
		   Safefree(nums);
		   nums = NULL;
	  }
	}

SV* matrix(...) 
CODE:
    /* Basic check: must have an even number of arguments for key => value */
    if (items % 2 != 0) {
        croak("Usage: matrix(data => [...], nrow => $n, ncol => $m, byrow => $bool)");
    }
    SV*restrict data_sv = NULL;
    UV  nrow = 0, ncol = 0;
    bool byrow = FALSE, nrow_set = FALSE, ncol_set = FALSE;
    /* Parse named arguments */
    for (I32 i = 0; i < items; i += 2) {
        char*restrict key = SvPV_nolen(ST(i));
        SV*restrict val   = ST(i + 1);
        if (strEQ(key, "data")) {
            data_sv = val;
        } else if (strEQ(key, "nrow")) {
            nrow = (UV)SvUV(val);
            nrow_set = TRUE;
        } else if (strEQ(key, "ncol")) {
            ncol = (UV)SvUV(val);
            ncol_set = TRUE;
        } else if (strEQ(key, "byrow")) {
            byrow = SvTRUE(val);
        } else {
            croak("Unknown option: %s", key);
        }
    }
    /* Validate data input */
    if (!data_sv || !SvROK(data_sv) || SvTYPE(SvRV(data_sv)) != SVt_PVAV) {
        croak("The 'data' option must be an array reference (e.g. data => [1..6])");
    }
    AV*restrict data_av = (AV*)SvRV(data_sv);
    UV  data_len = (UV)(av_top_index(data_av) + 1);
    if (data_len == 0) {
        croak("Data array cannot be empty");
    }
    /* R-style dimension inference */
    if (!nrow_set && !ncol_set) {
        nrow = data_len;
        ncol = 1;
    } else if (nrow_set && !ncol_set) {
        ncol = (data_len + nrow - 1) / nrow;
    } else if (!nrow_set && ncol_set) {
        nrow = (data_len + ncol - 1) / ncol;
    }
    /* Final safety check for dimensions */
    if (nrow == 0 || ncol == 0) {
        croak("Dimensions must be greater than 0");
    }
    /* Create the matrix (Array of Arrays) */
    AV*restrict result_av = newAV();
    av_extend(result_av, nrow - 1);
    UV r, c, i;/* Use unsigned types for counters to prevent negative indexing */
    AV**restrict row_ptrs = (AV**restrict)safemalloc(nrow * sizeof(AV*)); /* Pre-allocate row pointers */
    for (r = 0; r < nrow; r++) {
        row_ptrs[r] = newAV();
        av_extend(row_ptrs[r], ncol - 1);
        av_push(result_av, newRV_noinc((SV*)row_ptrs[r]));
    }
    /* Fill the matrix */
    UV total_cells = nrow * ncol;
    for (i = 0; i < total_cells; i++) {
        /* Vector recycling logic */
        SV**restrict fetched = av_fetch(data_av, i % data_len, 0);
        SV*restrict val = fetched ? newSVsv(*fetched) : newSV(0);
        if (byrow) {
            r = i / ncol;
            c = i % ncol;
        } else {
            r = i % nrow;
            c = i / nrow;
        }
        av_store(row_ptrs[r], c, val);
    }
    safefree(row_ptrs);
    RETVAL = newRV_noinc((SV*)result_av);
OUTPUT:
    RETVAL

SV* lm(...)
    CODE:
    {
        if (items % 2 != 0) 
            croak("Usage: lm(formula => 'mpg ~ wt * hp', data => \\%mtcars)");
        const char *restrict formula  = NULL;
        SV *restrict data_sv = NULL;
        for (I32 i = 0; i < items; i += 2) {
            const char *restrict key = SvPV_nolen(ST(i));
            SV *restrict val = ST(i + 1);
            if      (strEQ(key, "formula")) formula = SvPV_nolen(val);
            else if (strEQ(key, "data"))    data_sv = val;
            else croak("lm: unknown argument '%s'", key);
        }        
        if (!formula) croak("lm: formula is required");
        if (!data_sv || !SvROK(data_sv)) croak("lm: data is required and must be a reference");
        /* --- 1. Clean and Strip the Formula --- */
        char f_cpy[512];
        char *restrict src = (char*)formula;
        char *restrict dst = f_cpy;
        while (*src && (dst - f_cpy < 511)) {
            if (!isspace(*src)) { *dst++ = *src; }
            src++;
        }
        *dst = '\0';

        char *restrict tilde = strchr(f_cpy, '~');
        if (!tilde) croak("lm: invalid formula, missing '~'");
        *tilde = '\0';
        
        char *restrict lhs = f_cpy;         /* Response Variable */
        char *restrict rhs = tilde + 1;     /* Predictor Expression */

        /* --- 2. Tokenize Terms and Resolve Interaction Macros --- */
        char terms[128][128];
        unsigned int num_terms = 0;
        
        bool has_intercept = true;
        if (strstr(rhs, "-1")) has_intercept = false;
        
        if (has_intercept) {
            strcpy(terms[num_terms++], "Intercept");
        }
        /* Tokenize by the '+' operator */
        char *restrict chunk = strtok(rhs, "+");
        while (chunk != NULL) {
            if (strcmp(chunk, "1") == 0 || strcmp(chunk, "-1") == 0) {
                chunk = strtok(NULL, "+");
                continue;
            }
            char *restrict star = strchr(chunk, '*');
            if (star) {
                *star = '\0';
                char *restrict left = chunk;
                char *restrict right = star + 1;
                
                /* R formula rule: strip ^N unless wrapped in I() */
                char *restrict c_left = strchr(left, '^'); 
                if (c_left && strncmp(left, "I(", 2) != 0) *c_left = '\0';
                char *restrict c_right = strchr(right, '^'); 
                if (c_right && strncmp(right, "I(", 2) != 0) *c_right = '\0';
                
                strcpy(terms[num_terms++], left);
                strcpy(terms[num_terms++], right);
                sprintf(terms[num_terms++], "%s:%s", left, right);
            } else {
                char *restrict c_chunk = strchr(chunk, '^'); 
                if (c_chunk && strncmp(chunk, "I(", 2) != 0) *c_chunk = '\0';
                strcpy(terms[num_terms++], chunk);
            }
            chunk = strtok(NULL, "+");
        }

        char uniq_terms[128][128];
        unsigned int num_uniq = 0;
        for (unsigned int i = 0; i < num_terms; i++) {
            bool found = false;
            for (unsigned int j = 0; j < num_uniq; j++) {
                if (strcmp(terms[i], uniq_terms[j]) == 0) { found = true; break; }
            }
            if (!found) {
                strcpy(uniq_terms[num_uniq++], terms[i]);
            }
        }
        unsigned int p = num_uniq, n = 0;
        /* --- 3. Structure Row Iterators for Hash of Hashes --- */
        char **restrict row_names = NULL;
        HV **restrict row_hashes = NULL;
        HV *restrict data_hoa = NULL;

        SV*restrict ref = SvRV(data_sv);
        if (SvTYPE(ref) == SVt_PVHV) {
            HV*restrict hv = (HV*)ref;
            if (hv_iterinit(hv) == 0) croak("lm: Data hash is empty");
            HE*restrict entry = hv_iternext(hv);
            if (entry) {
                SV* val = hv_iterval(hv, entry);
                if (SvROK(val) && SvTYPE(SvRV(val)) == SVt_PVAV) {
                    /* Fallback: HoA mode */
                    data_hoa = hv;
                    n = av_len((AV*)SvRV(val)) + 1;
                    Newx(row_names, n, char*);
                    for(unsigned int i = 0; i < n; i++) {
                        char buf[32]; snprintf(buf, sizeof(buf), "%u", i+1);
                        row_names[i] = savepv(buf);
                    }
                } else if (SvROK(val) && SvTYPE(SvRV(val)) == SVt_PVHV) {
                    /* Primary: HoH mode (%mtcars keys are "Mazda RX4") */
                    n = hv_iterinit(hv);
                    Newx(row_names, n, char*);
                    Newx(row_hashes, n, HV*);
                    size_t i = 0;
                    while ((entry = hv_iternext(hv))) {
                        I32 len;
                        char *restrict key = hv_iterkey(entry, &len);
                        row_names[i] = savepv(key);
                        row_hashes[i] = (HV*)SvRV(hv_iterval(hv, entry));
                        i++;
                    }
                } else {
                    croak("lm: Hash values must be ArrayRefs (HoA) or HashRefs (HoH)");
                }
            }
        } else if (SvTYPE(ref) == SVt_PVAV) {
            /* Fallback: AoH mode */
            AV*restrict av = (AV*)ref;
            n = av_len(av) + 1;
            Newx(row_names, n, char*);
            Newx(row_hashes, n, HV*);
            for (size_t i = 0; i < n; i++) {
                SV**restrict val = av_fetch(av, i, 0);
                if (val && SvROK(*val) && SvTYPE(SvRV(*val)) == SVt_PVHV) {
                    row_hashes[i] = (HV*)SvRV(*val);
                    char buf[32]; snprintf(buf, sizeof(buf), "%lu", i + 1);
                    row_names[i] = savepv(buf);
                } else {
                    for (size_t x = 0; x < i; x++) Safefree(row_names[x]);
                    Safefree(row_names); Safefree(row_hashes);
                    croak("lm: Array values must be HashRefs (AoH)");
                }
            }
        } else {
            croak("lm: Data must be an Array or Hash reference");
        }

        double *restrict X; Newx(X, n * p, double);
        double *restrict Y; Newx(Y, n, double);

        for (size_t i = 0; i < n; i++) {
            Y[i] = evaluate_term(data_hoa, row_hashes, i, lhs);
            for (size_t j = 0; j < p; j++) {
                if (strcmp(uniq_terms[j], "Intercept") == 0) {
                    X[i * p + j] = 1.0;
                } else {
                    X[i * p + j] = evaluate_term(data_hoa, row_hashes, i, uniq_terms[j]);
                }
            }
        }
        /* --- 4. Ordinary Least Squares (OLS) Math --- */
        double *restrict XtX; Newxz(XtX, p * p, double);
        for (size_t i = 0; i < p; i++) {
            for (size_t j = 0; j < p; j++) {
                double sum = 0.0;
                for (size_t k = 0; k < n; k++) {
                    sum += X[k * p + i] * X[k * p + j];
                }
                XtX[i * p + j] = sum;
            }
        }
        double *restrict XtY; Newxz(XtY, p, double);
        for (size_t i = 0; i < p; i++) {
            double sum = 0.0;
            for (size_t k = 0; k < n; k++) {
                sum += X[k * p + i] * Y[k];
            }
            XtY[i] = sum;
        }

        if (invert_matrix_gj(XtX, p) != 0) {
            Safefree(X); Safefree(Y); Safefree(XtX); Safefree(XtY);
            if (row_hashes) Safefree(row_hashes);
            for (size_t i = 0; i < n; i++) Safefree(row_names[i]);
            Safefree(row_names);
            croak("lm: computation failed (singular/rank-deficient matrix)");
        }

        double *restrict beta; Newxz(beta, p, double);
        for (unsigned int i = 0; i < p; i++) {
            double sum = 0.0;
            for (unsigned int j = 0; j < p; j++) {
                sum += XtX[i * p + j] * XtY[j];
            }
            beta[i] = sum;
        }

        /* --- 5. Assemble Return Structure mirroring R's `lm` list --- */
        HV *restrict res_hv = newHV();
        HV *restrict coef_hv = newHV();
        HV *restrict fitted_hv = newHV();
        HV *restrict resid_hv = newHV();

        for (unsigned int j = 0; j < p; j++) {
            hv_store(coef_hv, uniq_terms[j], strlen(uniq_terms[j]), newSVnv(beta[j]), 0);
        }
        
        for (unsigned int i = 0; i < n; i++) {
            double y_hat = 0.0;
            for (unsigned int j = 0; j < p; j++) {
                y_hat += X[i * p + j] * beta[j];
            }
            double res = Y[i] - y_hat;
            
            /* Add outputs keyed by explicit rowname ("Mazda RX4" => ...) */
            hv_store(fitted_hv, row_names[i], (I32)strlen(row_names[i]), newSVnv(y_hat), 0);
            hv_store(resid_hv, row_names[i], (I32)strlen(row_names[i]), newSVnv(res), 0);
            Safefree(row_names[i]);
        }
        Safefree(row_names);

        hv_store(res_hv, "coefficients",  12, newRV_noinc((SV*)coef_hv), 0);
        hv_store(res_hv, "fitted.values", 13, newRV_noinc((SV*)fitted_hv), 0);
        hv_store(res_hv, "residuals",      9, newRV_noinc((SV*)resid_hv), 0);
        hv_store(res_hv, "df.residual",   11, newSVuv(n - p), 0);
        hv_store(res_hv, "rank",          4,  newSVuv(p), 0);

        /* Clean up */
        Safefree(X); 
        Safefree(Y); 
        Safefree(XtX); 
        Safefree(XtY); 
        Safefree(beta);
        if (row_hashes) Safefree(row_hashes);

        RETVAL = newRV_noinc((SV*)res_hv);
    }
    OUTPUT:
        RETVAL

SV* rnorm(...)
	CODE:
	{
	  if (items % 2 != 0)
		   croak("Usage: rnorm(n => 10, mean => 0, sd => 1) - must be even key/value pairs");

	  /* --- Parse named arguments from the flat stack --- */
	  ssize_t n = 0;
	  double mean = 0.0, sd = 1.0;
	  for (I32 i = 0; i < items; i += 2) {
		   const char* restrict key = SvPV_nolen(ST(i));
		   SV* restrict val = ST(i + 1);

		   if      (strEQ(key, "n"))      n    = (unsigned int)SvUV(val);
		   else if (strEQ(key, "mean"))   mean = SvNV(val);
		   else if (strEQ(key, "sd"))     sd   = SvNV(val);
		   else croak("rnorm: unknown argument '%s'", key);
	  }
	  if (sd < 0.0) croak("rnorm: standard deviation must be non-negative");
	  AV *restrict result_av = newAV();
	  if (n > 0) {
		   av_extend(result_av, n - 1);
		   /* Generate random normals using the Box-Muller transform */
		   for (size_t i = 0; i < n; ) {
		       double u, v, s;
		       do {
		           /* Drand01() hooks into Perl's internal PRNG, respecting Perl's srand() */
		           u = 2.0 * Drand01() - 1.0;
		           v = 2.0 * Drand01() - 1.0;
		           s = u * u + v * v;
		       } while (s >= 1.0 || s == 0.0);
		       double mul = sqrt(-2.0 * log(s) / s);
		       /* Box-Muller generates two independent values per iteration */
		       av_store(result_av, i++, newSVnv(mean + sd * u * mul));
		       if (i < n) {
		           av_store(result_av, i++, newSVnv(mean + sd * v * mul));
		       }
		   }
	  }
	  RETVAL = newRV_noinc((SV*)result_av);
	}
	OUTPUT:
	  RETVAL

PROTOTYPES: DISABLE

SV*
fisher_test(data_ref, conf_level = 0.95)
    SV* data_ref
    double conf_level
    PREINIT:
        size_t a = 0, b = 0, c = 0, d = 0;
        double p_val, mle_or, ci_low, ci_high;
        HV*restrict ret_hash;
        AV*restrict ci_array;
        HV*restrict est_hash;
    CODE:
        if (!SvROK(data_ref)) croak("fisher_test requires a reference to an Array or Hash");
        SV*restrict deref = SvRV(data_ref);
        /* Fast Path: 2D Array / AoA */
        if (SvTYPE(deref) == SVt_PVAV) {
            AV*restrict outer = (AV*)deref;
            if (av_len(outer) != 1) croak("Outer array must have exactly 2 rows");
            
            SV**restrict row1_ptr = av_fetch(outer, 0, 0);
            SV**restrict row2_ptr = av_fetch(outer, 1, 0);
            
            if (row1_ptr && row2_ptr && SvROK(*row1_ptr) && SvROK(*row2_ptr)) {
                AV*restrict row1 = (AV*)SvRV(*row1_ptr);
                AV*restrict row2 = (AV*)SvRV(*row2_ptr);
                
                a = SvIV(*av_fetch(row1, 0, 0));
                b = SvIV(*av_fetch(row1, 1, 0));
                c = SvIV(*av_fetch(row2, 0, 0));
                d = SvIV(*av_fetch(row2, 1, 0));
            } else {
                croak("Invalid 2D Array structure");
            }
        } 
        /* Complex Path: 2D Hash (Because Perl hash order is random, we extract the first 4 integers found in the nested hashes) */
        else if (SvTYPE(deref) == SVt_PVHV) {
            HV*restrict outer = (HV*)deref;
            HE*restrict outer_entry;
            I32 outer_len;
            size_t val_count = 0;
            int vals[4] = {0, 0, 0, 0};
            
            hv_iterinit(outer);
            while ((outer_entry = hv_iternext(outer))) {
                SV* inner_sv = hv_iterval(outer, outer_entry);
                if (SvROK(inner_sv) && SvTYPE(SvRV(inner_sv)) == SVt_PVHV) {
                    HV* inner = (HV*)SvRV(inner_sv);
                    HE* inner_entry;
                    hv_iterinit(inner);
                    while ((inner_entry = hv_iternext(inner)) && val_count < 4) {
                        vals[val_count++] = SvIV(hv_iterval(inner, inner_entry));
                    }
                }
            }
            if (val_count != 4) croak("2D Hash must contain exactly 4 values");
            a = vals[0]; b = vals[1]; c = vals[2]; d = vals[3];
        } else {
            croak("Input must be a 2D Array or 2D Hash");
        }
        /* Perform Calculations */
        p_val = exact_p_value(a, b, c, d);
        calculate_exact_stats(a, b, c, d, conf_level, &mle_or, &ci_low, &ci_high);

        /* Construct the Return HashRef purely in C */
        ret_hash = newHV();
        hv_stores(ret_hash, "p_value", newSVnv(p_val));
        hv_stores(ret_hash, "method", newSVpv("Fisher's Exact Test for Count Data", 0));
        
        ci_array = newAV();
        av_push(ci_array, newSVnv(ci_low));
        av_push(ci_array, newSVnv(ci_high));
        hv_stores(ret_hash, "conf_int", newRV_noinc((SV*)ci_array));
        
        est_hash = newHV();
        hv_stores(est_hash, "odds ratio", newSVnv(mle_or));
        hv_stores(ret_hash, "estimate", newRV_noinc((SV*)est_hash));
        /* Return the HashRef */
        RETVAL = newRV_noinc((SV*)ret_hash);
    OUTPUT:
        RETVAL
