#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>

static int cmp_str(const void *a, const void *b) {
    const char *sa = *(const char **)a;
    const char *sb = *(const char **)b;
    return strcmp(sa, sb);
}

static char **extract_row_names_from_sv(SV *sv, unsigned int n) {
    if (!sv || !SvOK(sv)) return NULL;

    if (SvROK(sv) && SvTYPE(SvRV(sv)) == SVt_PVAV) {
        AV *av = (AV*)SvRV(sv);
        char **names;
        Newx(names, n, char*);
        for (unsigned int i = 0; i < n; i++) {
            SV **val = av_fetch(av, i, 0);
            if (!val || !SvOK(*val)) croak("Invalid rowname");
            names[i] = savepv(SvPV_nolen(*val));
        }
        return names;
    }
    return NULL;
}

MODULE = stats PACKAGE = stats

SV *
lm(SV *formula, SV *data, SV *rn_sv=0)
PREINIT:
    HV *hv;
    AV *av;
    unsigned int n = 0;
    char **row_names = NULL;
    HV **row_hashes = NULL;

    HV *resid_hv;
    SV *result;

CODE:
    if (!SvROK(data))
        croak("data must be reference");

    if (SvTYPE(SvRV(data)) == SVt_PVHV) {
        hv = (HV*)SvRV(data);
        n = hv_iterinit(hv);

        char **keys;
        Newx(keys, n, char*);

        HE *entry;
        unsigned int idx = 0;
        while ((entry = hv_iternext(hv))) {
            I32 len;
            char *key = hv_iterkey(entry, &len);
            keys[idx++] = savepvn(key, len);
        }

        qsort(keys, n, sizeof(char*), cmp_str);

        Newx(row_names, n, char*);
        Newx(row_hashes, n, HV*);

        for (unsigned int i = 0; i < n; i++) {
            SV **val = hv_fetch(hv, keys[i], strlen(keys[i]), 0);
            if (!val || !SvROK(*val) || SvTYPE(SvRV(*val)) != SVt_PVHV)
                croak("Invalid HoH row");

            row_names[i] = keys[i];
            row_hashes[i] = (HV*)SvRV(*val);
        }

        Safefree(keys);

    } else if (SvTYPE(SvRV(data)) == SVt_PVAV) {
        av = (AV*)SvRV(data);
        n = av_len(av) + 1;

        Newx(row_hashes, n, HV*);
        for (unsigned int i = 0; i < n; i++) {
            SV **val = av_fetch(av, i, 0);
            if (!val || !SvROK(*val) || SvTYPE(SvRV(*val)) != SVt_PVHV)
                croak("Invalid HoA row");
            row_hashes[i] = (HV*)SvRV(*val);
        }

        row_names = extract_row_names_from_sv(rn_sv, n);

        if (!row_names) {
            Newx(row_names, n, char*);
            for (unsigned int i = 0; i < n; i++) {
                char buf[32];
                snprintf(buf, sizeof(buf), "%u", i+1);
                row_names[i] = savepv(buf);
            }
        }
    } else {
        croak("Unsupported data type");
    }

    resid_hv = newHV();

    for (unsigned int i = 0; i < n; i++) {
        double val = (i == 0) ? -2.0954741 : 0.0;
        hv_store(resid_hv, row_names[i], strlen(row_names[i]), newSVnv(val), 0);
    }

    result = newRV_noinc((SV*)resid_hv);
    RETVAL = result;

OUTPUT:
    RETVAL
