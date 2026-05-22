/* =========================================================================
 * oneway_test  –  Welch / classic one-way ANOVA
 *
 * ── Mode 1: hash of groups (original behaviour) ───────────────────────────
 *
 *   my $res = oneway_test(\%groups);
 *   my $res = oneway_test(\%groups, var_equal => 1);
 *
 *   \%groups  – keys are group labels, values are array refs of numbers.
 *               Every group must have >= 2 observations.
 *
 * ── Mode 2: formula – response ~ factor ───────────────────────────────────
 *
 *   my $res = oneway_test(\%data, formula => "yield ~ ctrl");
 *   my $res = oneway_test(\%data, formula => "yield ~ ctrl", var_equal => 1);
 *
 *   \%data must contain two keys matching the formula:
 *     "yield" => [ numeric response values ... ]
 *     "ctrl"  => [ group labels (strings or numbers, same length) ... ]
 *
 *   This mirrors R's:
 *     my_data <- stack(list(yield = yield, ctrl = ctrl))
 *     oneway.test(Value ~ Group, data = my_data)
 *
 *   Absence of a formula argument falls back to Mode 1 automatically.
 *
 * ── Return value (both modes) ─────────────────────────────────────────────
 *
 *   Hash ref with keys:
 *     statistic   => F value
 *     num_df      => numerator degrees of freedom   (k − 1)
 *     denom_df    => denominator degrees of freedom
 *     p_value     => upper-tail p-value  P(F ≥ statistic)
 *     method      => description string
 *     k           => number of groups
 *     n           => total observations
 *     formula     => "response ~ factor"  (only present in Mode 2)
 *
 * =========================================================================
 * Integration: drop the C block above "--- XS SECTION ---", and the XS
 * block inside your MODULE … PACKAGE … PREFIX = section.
 * =========================================================================
 */

/* -----------------------------------------------------------------------
 * C HELPERS  (place above "--- XS SECTION ---")
 * ----------------------------------------------------------------------- */

/* ── OneWayResult struct ─────────────────────────────────────────────── */
typedef struct {
    double  statistic;
    double  num_df;
    double  denom_df;
    double  p_value;
    int     k;          /* number of groups          */
    IV      n;          /* total observations        */
    int     var_equal;  /* 0 = Welch, 1 = classic    */
} OneWayResult;

/* ── c_oneway_test ───────────────────────────────────────────────────────
 *
 *  data      – flat C array of all observations, groups concatenated
 *  sizes     – n_i for each group (length k)
 *  k         – number of groups
 *  var_equal – 0 = Welch (default), 1 = classic equal-variance F-test
 *
 *  Mirrors R's oneway.test() arithmetic exactly.
 *  Calls pf() which must be declared elsewhere in the .xs file.
 * ----------------------------------------------------------------------- */
static OneWayResult
c_oneway_test(const double *restrict data,
              const size_t *restrict sizes,
              size_t k,
              int var_equal)
{
    OneWayResult res;
    res.var_equal = var_equal;
    res.k         = (int)k;

    double *restrict n_i = (double *)safemalloc(k * sizeof(double));
    double *restrict m_i = (double *)safemalloc(k * sizeof(double));
    double *restrict v_i = (double *)safemalloc(k * sizeof(double));

    size_t offset = 0;
    IV total_n = 0;
    for (size_t g = 0; g < k; g++) {
        size_t ng  = sizes[g];
        n_i[g]     = (double)ng;
        total_n   += (IV)ng;

        double sum = 0.0;
        for (size_t i = 0; i < ng; i++) sum += data[offset + i];
        double mean = sum / (double)ng;
        m_i[g] = mean;

        double ss = 0.0;
        for (size_t i = 0; i < ng; i++) {
            double d = data[offset + i] - mean;
            ss += d * d;
        }
        v_i[g] = ss / (double)(ng - 1);   /* ng >= 2 guaranteed by caller */
        offset += ng;
    }

    res.n = total_n;

    /* grand mean (simple average over all obs; used only by classic branch) */
    double grand_mean = 0.0;
    for (IV i = 0; i < (IV)total_n; i++) grand_mean += data[i];
    grand_mean /= (double)total_n;

    double df1 = (double)(k - 1);

    if (var_equal) {
        /* ── Classic one-way ANOVA ─────────────────────────────────────── *
         *  F = [Σ n_i·(m_i − ȳ)² / (k−1)]  /  [Σ (n_i−1)·v_i / (n−k)] *
         * ─────────────────────────────────────────────────────────────── */
        double ssbg = 0.0, sswg = 0.0;
        for (size_t g = 0; g < k; g++) {
            double dm = m_i[g] - grand_mean;
            ssbg += n_i[g] * dm * dm;
            sswg += (n_i[g] - 1.0) * v_i[g];
        }
        double df2   = (double)(total_n - (IV)k);
        res.statistic = (ssbg / df1) / (sswg / df2);
        res.num_df    = df1;
        res.denom_df  = df2;
    } else {
        /* ── Welch one-way (heteroscedastic) ───────────────────────────── *
         *  w_i  = n_i / v_i                                               *
         *  W    = Σ w_i                                                   *
         *  m̃    = Σ(w_i·m_i) / W          (weighted grand mean)           *
         *  tmp  = Σ[(1 − w_i/W)² / (n_i−1)] / (k²−1)                    *
         *  F    = Σ[w_i·(m_i − m̃)²] / [(k−1)·(1 + 2·(k−2)·tmp)]        *
         *  df2  = 1 / (3·tmp)                                             *
         * ─────────────────────────────────────────────────────────────── */
        double *restrict w_i = (double *)safemalloc(k * sizeof(double));

        double sum_w = 0.0;
        for (size_t g = 0; g < k; g++) { w_i[g] = n_i[g] / v_i[g]; sum_w += w_i[g]; }

        double wgrand = 0.0;
        for (size_t g = 0; g < k; g++) wgrand += w_i[g] * m_i[g];
        wgrand /= sum_w;

        double tmp = 0.0;
        for (size_t g = 0; g < k; g++) {
            double t = 1.0 - w_i[g] / sum_w;
            tmp += (t * t) / (n_i[g] - 1.0);
        }
        tmp /= ((double)k * (double)k - 1.0);   /* k² − 1 */

        double num = 0.0;
        for (size_t g = 0; g < k; g++) {
            double dm = m_i[g] - wgrand;
            num += w_i[g] * dm * dm;
        }

        res.statistic = num / (df1 * (1.0 + 2.0 * (double)(k - 2) * tmp));
        res.num_df    = df1;
        res.denom_df  = (tmp > 0.0) ? (1.0 / (3.0 * tmp)) : 1e300;

        Safefree(w_i);
    }

    /* upper-tail p-value  P(F ≥ statistic) */
    res.p_value = pf(res.statistic, res.num_df, res.denom_df, 0 /*lower*/, 0 /*log*/);

    Safefree(n_i);
    Safefree(m_i);
    Safefree(v_i);
    return res;
}

/* ── parse_formula ───────────────────────────────────────────────────────
 *
 *  Splits "response ~ factor" into two NUL-terminated, heap-allocated
 *  strings.  Leading/trailing whitespace is stripped from each side.
 *  Returns 1 on success, 0 on failure (malformed / missing '~').
 *  Caller must Safefree() both *lhs and *rhs on success.
 * ----------------------------------------------------------------------- */
static int
parse_formula(const char *formula, char **lhs, char **rhs)
{
    const char *tilde = strchr(formula, '~');
    if (!tilde) return 0;

    /* left-hand side: trim trailing whitespace */
    const char *l_start = formula;
    const char *l_end   = tilde - 1;
    while (l_end >= l_start && isspace((unsigned char)*l_end)) l_end--;
    if (l_end < l_start) return 0;             /* empty LHS */

    /* right-hand side: trim leading whitespace */
    const char *r_start = tilde + 1;
    while (*r_start && isspace((unsigned char)*r_start)) r_start++;
    const char *r_end = r_start + strlen(r_start) - 1;
    while (r_end >= r_start && isspace((unsigned char)*r_end)) r_end--;
    if (r_end < r_start) return 0;             /* empty RHS */

    size_t llen = (size_t)(l_end - l_start + 1);
    size_t rlen = (size_t)(r_end - r_start + 1);

    *lhs = (char *)safemalloc(llen + 1);
    *rhs = (char *)safemalloc(rlen + 1);
    memcpy(*lhs, l_start, llen); (*lhs)[llen] = '\0';
    memcpy(*rhs, r_start, rlen); (*rhs)[rlen] = '\0';
    return 1;
}

/* ── build_groups_from_formula ───────────────────────────────────────────
 *
 *  Takes parallel response[] and label[] arrays (each length n) and
 *  partitions them into groups, filling:
 *    out_flat[]  – observations sorted into contiguous group blocks
 *    out_sizes[] – number of observations per group  (caller allocates n
 *                  slots for both; actual group count returned via *out_k)
 *
 *  Group identity is the string representation of each label element
 *  (SvPV_nolen), so integer 0 and string "0" are the same group.
 *  Groups are ordered by first appearance in label[], matching R's
 *  factor level ordering from stack().
 *
 *  Returns 1 on success; 0 if any validation error (sets errbuf).
 * ----------------------------------------------------------------------- */
#define OWT_MAX_GROUPS 1024   /* sane ceiling; ANOVA with >1024 groups is absurd */

static int
build_groups_from_formula(pTHX_
                          AV *restrict response_av,
                          AV *restrict label_av,
                          double *restrict out_flat,
                          size_t *restrict out_sizes,
                          size_t *out_k,
                          char *errbuf,
                          size_t errbuf_len)
{
    IV n = av_len(response_av) + 1;
    IV nl = av_len(label_av)   + 1;

    if (n != nl) {
        snprintf(errbuf, errbuf_len,
            "formula: response length (%"IVdf") != factor length (%"IVdf")",
            n, nl);
        return 0;
    }
    if (n < 2) {
        snprintf(errbuf, errbuf_len, "formula: need at least 2 observations");
        return 0;
    }

    /* ── discover unique group labels in order of first appearance ─── */
    /* We store pointers into a heap-allocated label string table.       */
    char  **group_names  = (char **)safemalloc(OWT_MAX_GROUPS * sizeof(char *));
    size_t  ngroups      = 0;
    IV     *obs_group    = (IV *)safemalloc((size_t)n * sizeof(IV));
                                             /* maps obs index → group index */

    for (IV i = 0; i < n; i++) {
        SV **lsv = av_fetch(label_av, i, 0);
        const char *label = (lsv && *lsv) ? SvPV_nolen(*lsv) : "";

        /* linear scan for existing group (k is small, O(n·k) is fine) */
        IV gidx = -1;
        for (size_t g = 0; g < ngroups; g++) {
            if (strEQ(group_names[g], label)) { gidx = (IV)g; break; }
        }
        if (gidx < 0) {
            if (ngroups >= OWT_MAX_GROUPS) {
                snprintf(errbuf, errbuf_len,
                    "formula: too many distinct groups (max %d)", OWT_MAX_GROUPS);
                Safefree(group_names);
                Safefree(obs_group);
                return 0;
            }
            /* new group: copy the label string */
            size_t lablen = strlen(label);
            group_names[ngroups] = (char *)safemalloc(lablen + 1);
            memcpy(group_names[ngroups], label, lablen + 1);
            gidx = (IV)ngroups++;
        }
        obs_group[i] = gidx;
    }

    if (ngroups < 2) {
        snprintf(errbuf, errbuf_len,
            "formula: need at least 2 distinct groups, found %zu", ngroups);
        for (size_t g = 0; g < ngroups; g++) Safefree(group_names[g]);
        Safefree(group_names);
        Safefree(obs_group);
        return 0;
    }

    /* ── count per-group sizes ─────────────────────────────────────── */
    memset(out_sizes, 0, ngroups * sizeof(size_t));
    for (IV i = 0; i < n; i++) out_sizes[obs_group[i]]++;

    /* validate: every group needs >= 2 observations */
    for (size_t g = 0; g < ngroups; g++) {
        if (out_sizes[g] < 2) {
            snprintf(errbuf, errbuf_len,
                "formula: group '%s' has only %zu observation(s); need >= 2",
                group_names[g], out_sizes[g]);
            for (size_t gg = 0; gg < ngroups; gg++) Safefree(group_names[gg]);
            Safefree(group_names);
            Safefree(obs_group);
            return 0;
        }
    }

    /* ── fill flat output array in group order ─────────────────────── *
     *  We compute a running write-offset per group, then scatter.      *
     * ─────────────────────────────────────────────────────────────── */
    size_t *restrict write_pos = (size_t *)safemalloc(ngroups * sizeof(size_t));
    write_pos[0] = 0;
    for (size_t g = 1; g < ngroups; g++)
        write_pos[g] = write_pos[g - 1] + out_sizes[g - 1];

    for (IV i = 0; i < n; i++) {
        SV **rsv = av_fetch(response_av, i, 0);
        double val = (rsv && *rsv) ? SvNV(*rsv) : 0.0;
        size_t g   = (size_t)obs_group[i];
        out_flat[write_pos[g]++] = val;
    }

    *out_k = ngroups;

    /* ── clean up ──────────────────────────────────────────────────── */
    Safefree(write_pos);
    Safefree(obs_group);
    for (size_t g = 0; g < ngroups; g++) Safefree(group_names[g]);
    Safefree(group_names);
    return 1;
}
#undef OWT_MAX_GROUPS


/* -----------------------------------------------------------------------
 * XS WRAPPER  (place inside your MODULE … PACKAGE … PREFIX = block)
 *
 * Signatures:
 *   oneway_test(\%groups)
 *   oneway_test(\%groups, var_equal  => 1)
 *   oneway_test(\%data,   formula    => "response ~ factor")
 *   oneway_test(\%data,   formula    => "response ~ factor", var_equal => 1)
 * ----------------------------------------------------------------------- */

SV *
oneway_test(hashref, ...)
    SV *hashref
  PREINIT:
    HV          *in_hv;
    HE          *he;
    int          var_equal   = 0;
    const char  *formula_str = NULL;   /* NULL = Mode 1 (groups hash) */
    char        *lhs = NULL, *rhs = NULL;
    double      *flat  = NULL;
    size_t      *sizes = NULL;
    size_t       k     = 0;
    IV           total_n = 0;
    OneWayResult res;
    HV          *ret_hv;
    char         errbuf[512];

  CODE:
    /* ── parse named arguments ───────────────────────────────────────── */
    for (I32 ai = 1; ai + 1 < items; ai += 2) {
        const char *key = SvPV_nolen(ST(ai));
        SV         *val = ST(ai + 1);
        if (strEQ(key, "var_equal"))
            var_equal = SvTRUE(val) ? 1 : 0;
        else if (strEQ(key, "formula"))
            formula_str = SvPV_nolen(val);
    }

    /* ── validate hashref ────────────────────────────────────────────── */
    if (!SvROK(hashref) || SvTYPE(SvRV(hashref)) != SVt_PVHV)
        croak("oneway_test: first argument must be a hash reference");
    in_hv = (HV *)SvRV(hashref);

    if (formula_str != NULL) {
        /* ══════════════════════════════════════════════════════════════
         * MODE 2 – formula  "response ~ factor"
         *
         *  Parse the formula, look up both arrays in \%data,
         *  then delegate to build_groups_from_formula().
         * ══════════════════════════════════════════════════════════════ */
        if (!parse_formula(formula_str, &lhs, &rhs))
            croak("oneway_test: cannot parse formula '%s' — "
                  "expected 'response ~ factor'", formula_str);

        /* look up response array */
        SV **resp_svp = hv_fetch(in_hv, lhs, (I32)strlen(lhs), 0);
        if (!resp_svp || !*resp_svp || !SvROK(*resp_svp)
            || SvTYPE(SvRV(*resp_svp)) != SVt_PVAV)
            croak("oneway_test: formula LHS '%s' not found as an array ref "
                  "in the hash", lhs);

        /* look up factor/grouping array */
        SV **fact_svp = hv_fetch(in_hv, rhs, (I32)strlen(rhs), 0);
        if (!fact_svp || !*fact_svp || !SvROK(*fact_svp)
            || SvTYPE(SvRV(*fact_svp)) != SVt_PVAV)
            croak("oneway_test: formula RHS '%s' not found as an array ref "
                  "in the hash", rhs);

        AV *resp_av  = (AV *)SvRV(*resp_svp);
        AV *label_av = (AV *)SvRV(*fact_svp);
        IV  n        = av_len(resp_av) + 1;

        /* allocate worst-case buffers (at most n groups) */
        flat  = (double *)safemalloc((size_t)n * sizeof(double));
        sizes = (size_t *)safemalloc((size_t)n * sizeof(size_t));

        if (!build_groups_from_formula(aTHX_ resp_av, label_av,
                                       flat, sizes, &k,
                                       errbuf, sizeof errbuf)) {
            Safefree(flat);
            Safefree(sizes);
            Safefree(lhs);
            Safefree(rhs);
            croak("oneway_test: %s", errbuf);
        }

        for (size_t g = 0; g < k; g++) total_n += (IV)sizes[g];

    } else {
        /* ══════════════════════════════════════════════════════════════
         * MODE 1 – hash of groups  { label => \@observations, … }
         *
         *  Original behaviour, unchanged.
         * ══════════════════════════════════════════════════════════════ */
        k = (size_t)hv_iterinit(in_hv);
        if (k < 2)
            croak("oneway_test: need at least 2 groups, got %zu", k);

        /* first pass: sizes + total */
        sizes = (size_t *)safemalloc(k * sizeof(size_t));
        {
            size_t g = 0;
            while ((he = hv_iternext(in_hv)) != NULL) {
                SV *val = HeVAL(he);
                if (!SvROK(val) || SvTYPE(SvRV(val)) != SVt_PVAV)
                    croak("oneway_test: value for group '%s' is not an array ref",
                          HePV(he, PL_na));
                IV len = av_len((AV *)SvRV(val)) + 1;
                if (len < 2)
                    croak("oneway_test: group '%s' has fewer than 2 observations",
                          HePV(he, PL_na));
                sizes[g++] = (size_t)len;
                total_n   += (IV)len;
            }
        }

        /* second pass: fill flat */
        flat = (double *)safemalloc((size_t)total_n * sizeof(double));
        {
            size_t offset = 0;
            hv_iterinit(in_hv);
            while ((he = hv_iternext(in_hv)) != NULL) {
                AV *av  = (AV *)SvRV(HeVAL(he));
                IV  len = av_len(av) + 1;
                for (IV i = 0; i < len; i++) {
                    SV **svp = av_fetch(av, i, 0);
                    flat[offset++] = (svp && *svp) ? SvNV(*svp) : 0.0;
                }
            }
        }
    }

    /* ── run the arithmetic ──────────────────────────────────────────── */
    res = c_oneway_test(flat, sizes, k, var_equal);

    Safefree(flat);
    Safefree(sizes);
    if (lhs) Safefree(lhs);
    if (rhs) Safefree(rhs);

    /* ── build return hash ref ───────────────────────────────────────── */
    ret_hv = (HV *)sv_2mortal((SV *)newHV());

    hv_stores(ret_hv, "statistic", newSVnv(res.statistic));
    hv_stores(ret_hv, "num_df",    newSVnv(res.num_df));
    hv_stores(ret_hv, "denom_df",  newSVnv(res.denom_df));
    hv_stores(ret_hv, "p_value",   newSVnv(res.p_value));
    hv_stores(ret_hv, "k",         newSViv((IV)res.k));
    hv_stores(ret_hv, "n",         newSViv(res.n));
    hv_stores(ret_hv, "method",
        newSVpv(var_equal
            ? "One-way analysis of means"
            : "One-way analysis of means (not assuming equal variances)",
            0));

    /* only include 'formula' key when Mode 2 was used */
    if (formula_str)
        hv_stores(ret_hv, "formula", newSVpv(formula_str, 0));

    RETVAL = newRV((SV *)ret_hv);

  OUTPUT:
    RETVAL

/* -----------------------------------------------------------------------
 * END oneway_test
 * ----------------------------------------------------------------------- */
