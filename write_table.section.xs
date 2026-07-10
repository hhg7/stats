#ifndef _GNU_SOURCE
#define _GNU_SOURCE        // glibc / Linux
#endif
#ifndef __EXTENSIONS__
#define __EXTENSIONS__ 1   // Solaris/illumos: expose off64_t, sigjmp_buf under -std=c99
#endif
#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
#include "ppport.h"
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <strings.h>
#include <stdint.h> // uint64_t — harmless if perl.h already pulled it in
/*
XS words:
SvROK = scalar value reference is OK
*/
/* sample(): private splitmix64 PRNG

static bool contains_nondigit(pTHX_ SV *restrict sv) {
	if (!sv || !SvOK(sv)) return 0;
	STRLEN len;
	char *restrict s = SvPVbyte(sv, len);
	for (size_t i = 0; i < len; i++) {
	  if (!isdigit(s[i])) return 1;
	}
	return 0;
}

/* ---------------------------------------------------------------------------
 * print_string_row: emit one record.
 *
 * Quoting contract (matches the behaviors your tests pin down):
 *	- A field is quoted IFF it contains the separator string, a double
 *	  quote, or a newline / carriage return. Quoting is per-field, so in a
 *	  TSV "hello,world" stays bare while "tab\tin" becomes "tab	in".
 *	- Inside a quoted field, embedded double quotes are doubled:
 *	  p"q -> "p""q"   (RFC 4180 style)
 *	- A NULL or zero-length field prints as NOTHING between separators:
 *	  a,,c  -- never '' or "". A zero-length field cannot contain a
 *	  separator, quote, or newline, so it never needs quoting.
 *	- The separator is treated as a string (strstr), so multi-character
 *	  separators work; an empty separator never triggers quoting.
 *
 * Returns nothing; I/O errors surface on PerlIO_close at the call site.
*/
static void print_string_row(pTHX_ PerlIO *restrict fh,
	const char **restrict fields, size_t n, const char *restrict sep)
{
	const size_t sep_len = sep ? strlen(sep) : 0;
	for (size_t i = 0; i < n; i++) {
		if (i && sep_len) PerlIO_write(fh, sep, sep_len);
		const char *restrict f = fields[i];
		if (!f || !*f) continue; /* undef/empty -> print nothing */
		/* Does this field need quoting? */
		bool need_quotes = 0;
		if (strchr(f, '"') || strchr(f, '\n') || strchr(f, '\r')) {
			need_quotes = 1;
		} else if (sep_len && strstr(f, sep)) {
			need_quotes = 1;
		}
		if (!need_quotes) {
			PerlIO_write(fh, f, strlen(f));
		} else {
			PerlIO_putc(fh, '"');
			for (const char *restrict p = f; *p; p++) {
				if (*p == '"') PerlIO_putc(fh, '"'); /* double it */
				PerlIO_putc(fh, *p);
			}
			PerlIO_putc(fh, '"');
		}
	}
	PerlIO_putc(fh, '\n');
}

// Calculates the Regularized Upper Incomplete Gamma Function Q(a, x)
// This perfectly replicates R's pchisq(..., lower.tail=FALSE)
typedef int (*cs_cmp_fn)(pTHX_ void *restrict ctx, size_t i, size_t j);

/* Sort by a named column: pre-fetched cell SVs plus a numeric/string flag. */
typedef struct {
	SV **restrict vals;	/* borrowed cell SV* per row (NULL == missing) */
	unsigned short numeric;	/* 1 => compare with SvNV, 0 => compare with sv_cmp */
} cs_col_ctx;

/* Sort by a user comparator: per-row refs handed to $a/$b before each call. */
typedef struct {
	SV **restrict rows;	/* row ref per index (RV to HV) */
	CV  *restrict cv;	/* the comparator */
	SV  *a_sv;		/* scalar currently aliased to package $a */
	SV  *b_sv;		/* scalar currently aliased to package $b */
} cs_code_ctx;

MODULE = Stats::LikeR  PACKAGE = Stats::LikeR

void write_table(...)
PPCODE:
{
	SV *restrict data_sv = NULL;
	SV *restrict file_sv = NULL;
	unsigned int arg_idx = 0;
	// Mimic the Perl shift logic
	if (arg_idx < items && SvROK(ST(arg_idx))) {
		int type = SvTYPE(SvRV(ST(arg_idx)));
		if (type == SVt_PVHV || type == SVt_PVAV) {
			data_sv = ST(arg_idx);
			arg_idx++;
		}
	}
// Only consume a positional file argument if it is a plain string that is
// NOT one of the named option keys. Otherwise write_table(data=>..., file=>...)
// would grab the literal string "data" as the filename.
	if (arg_idx < items) {
		SV *restrict cand = ST(arg_idx);
		if (SvOK(cand) && !SvROK(cand)) {
			const char *restrict k = SvPV_nolen(cand);
			if (!(strEQ(k, "data") || strEQ(k, "file") || strEQ(k, "col.names") ||
				  strEQ(k, "row.names") || strEQ(k, "sep") || strEQ(k, "delim") ||
				  strEQ(k, "undef.val"))) {
				file_sv = cand;
				arg_idx++;
			}
		}
	}
	const char *restrict sep = ",";
	bool explicit_sep = 0; // Track if delimiter was manually specified
	// CHANGED: default undef cells to a true empty value ("") instead of NULL.
	// With print_string_row emitting zero-length fields bare (no quotes), an
	// undef cell now prints as nothing at all: a,,c -- not a,'',c or a,"",c.
	// 'undef.val' => 'NA' (etc.) still overrides this.
	const char *restrict undef_val = "";
	SV *restrict row_names_sv = sv_2mortal(newSViv(1));
	SV *restrict col_names_sv = NULL;
	// Read the remaining Hash-style arguments
	for (; arg_idx < items; arg_idx += 2) {
		if (arg_idx + 1 >= items) croak("write_table: Odd number of arguments passed");
		const char *restrict key = SvPV_nolen(ST(arg_idx));
		SV *restrict val = ST(arg_idx + 1);
		if (strEQ(key, "data")) data_sv = val;
		else if (strEQ(key, "col.names")) col_names_sv = val;
		else if (strEQ(key, "file")) file_sv = val;
		else if (strEQ(key, "row.names")) row_names_sv = val;
		// Check for either "sep" or "delim" and mark as explicitly provided
		else if (strEQ(key, "sep") || strEQ(key, "delim")) {
			sep = SvPV_nolen(val);
			explicit_sep = 1;
		}
		// FIX: 'undef.val' => undef used to call SvPV_nolen(&PL_sv_undef)
		// (warning + empty string by accident); make it explicit.
		else if (strEQ(key, "undef.val")) undef_val = SvOK(val) ? SvPV_nolen(val) : "";
		else croak("write_table: Unknown arguments passed: %s", key);
	}
	if (!data_sv || !SvROK(data_sv)) {
		croak("write_table: 'data' must be a HASH or ARRAY reference\n");
	}
	SV *restrict data_ref = SvRV(data_sv);
	if (SvTYPE(data_ref) != SVt_PVHV && SvTYPE(data_ref) != SVt_PVAV) {
		croak("write_table: 'data' must be a HASH or ARRAY reference\n");
	}
	if (!file_sv || !SvOK(file_sv)) croak("write_table: file name missing\n");
	const char *restrict file = SvPV_nolen(file_sv);
	// Auto-detect separator from file extension if not overridden
	if (!explicit_sep) {
		size_t file_len = strlen(file);
		if (file_len >= 4) {
			const char *restrict ext = file + file_len - 4;
			if (strEQ(ext, ".tsv") || strEQ(ext, ".TSV")) {
				sep = "\t";
			} else if (strEQ(ext, ".csv") || strEQ(ext, ".CSV")) {
				sep = ",";
			}
		}
	}
	if (col_names_sv && SvOK(col_names_sv)) {
		if (!SvROK(col_names_sv) || SvTYPE(SvRV(col_names_sv)) != SVt_PVAV) {
			croak("write_table: 'col.names' must be an ARRAY reference\n");
		}
	}
	bool is_hoh = 0, is_hoa = 0, is_aoh = 0, is_flat_hash = 0;
	AV *restrict rows_av = NULL;
// Validate Input Structures & Homogeneity
	if (SvTYPE(data_ref) == SVt_PVHV) {
		HV *restrict hv = (HV*)data_ref;
		if (hv_iterinit(hv) == 0) XSRETURN_EMPTY;
		HE *restrict entry = hv_iternext(hv);
		SV *restrict first_val = hv_iterval(hv, entry);

		if (!first_val) {
			croak("write_table: Invalid hash entry\n");
		}
// Check if top level values are scalars (Flat Hash)
		if (!SvROK(first_val)) {
			is_flat_hash = 1;
		} else {
			int first_type = SvTYPE(SvRV(first_val));
			if (first_type != SVt_PVHV && first_type != SVt_PVAV) {
				croak("write_table: Data values must be either all HASHes, all ARRAYs, or all scalars\n");
			}
			is_hoh = (first_type == SVt_PVHV);
			is_hoa = (first_type == SVt_PVAV);
		}
		hv_iterinit(hv);
		while ((entry = hv_iternext(hv))) {
			SV *restrict val = hv_iterval(hv, entry);
			if (is_flat_hash) {
				if (val && SvROK(val)) {
					croak("write_table: Mixed data types detected. Ensure all values are scalars for a flat hash.\n");
				}
			} else {
				if (!val || !SvROK(val) || SvTYPE(SvRV(val)) != (is_hoh ? SVt_PVHV : SVt_PVAV)) {
					croak("write_table: Mixed data types detected. Ensure all values are %s references.\n", is_hoh ? "HASH" : "ARRAY");
				}
			}
		}
		if (is_hoh) { // Rows are only explicitly pre-gathered for HOH
			rows_av = newAV();
			hv_iterinit(hv);
			while ((entry = hv_iternext(hv))) {
				av_push(rows_av, newSVsv(hv_iterkeysv(entry)));
			}
		}
	} else {
		AV *restrict av = (AV*)data_ref;
		if (av_len(av) < 0) XSRETURN_EMPTY;
		SV **restrict first_ptr = av_fetch(av, 0, 0);
		if (!first_ptr || !*first_ptr || !SvROK(*first_ptr) || SvTYPE(SvRV(*first_ptr)) != SVt_PVHV) {
			if (first_ptr && *first_ptr && SvROK(*first_ptr))
				croak("write_table: For ARRAY data, every element must be a HASH reference "
					  "(Array of Hashes); element 0 is a reference of type '%s'\n",
					  sv_reftype(SvRV(*first_ptr), 0));
			else if (first_ptr && *first_ptr && SvOK(*first_ptr))
				croak("write_table: For ARRAY data, every element must be a HASH reference "
					  "(Array of Hashes); element 0 is a non-reference scalar (value: '%s')\n",
					  SvPV_nolen(*first_ptr));
			else
				croak("write_table: For ARRAY data, every element must be a HASH reference "
					  "(Array of Hashes); element 0 is undef\n");
		}
// FIX: i was size_t while av_len() returns SSize_t; keep both signed.
		for (SSize_t i = 0; i <= av_len(av); i++) {
			SV **restrict ptr = av_fetch(av, i, 0);
			if (!ptr || !*ptr || !SvROK(*ptr) || SvTYPE(SvRV(*ptr)) != SVt_PVHV) {
				croak("write_table: Mixed data types detected in Array of Hashes. All elements must be HASH references.\n");
			}
		}
		is_aoh = 1;
	}
	PerlIO *restrict fh = PerlIO_open(file, "w");
	if (!fh) {
		// FIX: rows_av was leaked here when the open failed on HoH input.
		if (rows_av) SvREFCNT_dec(rows_av);
		croak("write_table: Could not open '%s' for writing", file);
	}
	AV *restrict headers_av = newAV();
	bool inc_rownames = (row_names_sv && SvTRUE(row_names_sv)) ? 1 : 0;
	const char *restrict rownames_col = NULL;
	// ----- Hash of Hashes -----
	if (is_hoh) {
		if (col_names_sv && SvOK(col_names_sv)) {
			AV *restrict c_av = (AV*)SvRV(col_names_sv);
			// FIX: i was size_t; av_len() == -1 on an empty col.names array
			// converted to SIZE_MAX and looped (effectively) forever.
			for (SSize_t i = 0; i <= av_len(c_av); i++) {
				SV **restrict c = av_fetch(c_av, i, 0);
				if (c && SvOK(*c)) av_push(headers_av, newSVsv(*c));
			}
		} else {
			HV *restrict col_map = newHV();
			hv_iterinit((HV*)data_ref);
			HE *restrict entry;
			while ((entry = hv_iternext((HV*)data_ref))) {
				HV *restrict inner = (HV*)SvRV(hv_iterval((HV*)data_ref, entry));
				hv_iterinit(inner);
				HE *restrict inner_entry;
				while ((inner_entry = hv_iternext(inner))) {
					hv_store_ent(col_map, hv_iterkeysv(inner_entry), newSViv(1), 0);
				}
			}
			unsigned num_cols = hv_iterinit(col_map);
			// FIX (UTF-8 safety): keep the key SVs (flags intact) and sort
			// them with sv_cmp instead of round-tripping through char*.
			for (unsigned i = 0; i < num_cols; i++) {
				HE *restrict ce = hv_iternext(col_map);
				av_push(headers_av, newSVsv(hv_iterkeysv(ce)));
			}
			if (num_cols > 1)
				sortsv(AvARRAY(headers_av), num_cols, Perl_sv_cmp);
			SvREFCNT_dec(col_map);
		}
		size_t num_headers = (size_t)(av_len(headers_av) + 1);
		const char **restrict header_row = safemalloc((num_headers + 1) * sizeof(char*));
		size_t h_idx = 0;
		if (inc_rownames) header_row[h_idx++] = "";
		// FIX: loop index was 'unsigned short int' -- silently wraps (and
		// loops forever) past 65535 columns. Use size_t like everywhere else.
		for (size_t i = 0; i < num_headers; i++) {
			SV **restrict h_ptr = av_fetch(headers_av, (SSize_t)i, 0);
			header_row[h_idx++] = (h_ptr && SvOK(*h_ptr)) ? SvPV_nolen(*h_ptr) : "";
		}
		print_string_row(aTHX_ fh, header_row, h_idx, sep);
		safefree(header_row);
		size_t num_rows = (size_t)(av_len(rows_av) + 1);
		// FIX (UTF-8/NUL safety): sort the key SVs themselves and look rows
		// up by SV (hv_fetch_ent) so UTF-8-flagged or NUL-containing outer
		// keys still match. sortsv+sv_cmp is plain string order, as before.
		sortsv(AvARRAY(rows_av), num_rows, Perl_sv_cmp);
		HV *restrict data_hv = (HV*)data_ref;
		const char **restrict row_data = safemalloc((num_headers + 1) * sizeof(char*));
		for (size_t i = 0; i < num_rows; i++) {
			size_t d_idx = 0;
			SV *restrict row_key_sv = *av_fetch(rows_av, (SSize_t)i, 0);
			if (inc_rownames) row_data[d_idx++] = SvPV_nolen(row_key_sv);
			HE *restrict inner_he = hv_fetch_ent(data_hv, row_key_sv, 0, 0);
			SV *restrict inner_sv = inner_he ? HeVAL(inner_he) : NULL;
			HV *restrict inner_hv = (inner_sv && SvROK(inner_sv)) ? (HV*)SvRV(inner_sv) : NULL;
			for (size_t j = 0; j < num_headers; j++) {
				SV **restrict h_ptr = av_fetch(headers_av, (SSize_t)j, 0);
				SV *restrict h_sv = (h_ptr && SvOK(*h_ptr)) ? *h_ptr : NULL;
				// FIX (UTF-8/NUL safety): fetch by SV, not by raw bytes
				HE *restrict cell_he = (inner_hv && h_sv) ? hv_fetch_ent(inner_hv, h_sv, 0, 0) : NULL;
				SV *restrict cell_sv = cell_he ? HeVAL(cell_he) : NULL;
				if (cell_sv && SvOK(cell_sv)) {
					if (SvROK(cell_sv)) {
						PerlIO_close(fh);
						safefree(row_data);
						if (headers_av) SvREFCNT_dec(headers_av);
						if (rows_av) SvREFCNT_dec(rows_av);
						croak("write_table: Cannot write nested reference types to table\n");
					}
					row_data[d_idx++] = SvPV_nolen(cell_sv);
				} else {
					row_data[d_idx++] = undef_val;
				}
			}
			print_string_row(aTHX_ fh, row_data, d_idx, sep);
		}
		safefree(row_data);
	// ----- Flat Hash -----
	} else if (is_flat_hash) {
		HV *restrict data_hv = (HV*)data_ref;
		if (col_names_sv && SvOK(col_names_sv)) {
			AV *restrict c_av = (AV*)SvRV(col_names_sv);
			for (SSize_t i = 0; i <= av_len(c_av); i++) {
				SV **restrict c = av_fetch(c_av, i, 0);
				if (c && SvOK(*c)) av_push(headers_av, newSVsv(*c));
			}
		} else {
			// FIX (UTF-8 safety): keep the key SVs (flags intact) and sort
			// them with sv_cmp instead of round-tripping through char*.
			unsigned int num_cols = hv_iterinit(data_hv);
			for (unsigned int i = 0; i < num_cols; i++) {
				HE *restrict ce = hv_iternext(data_hv);
				av_push(headers_av, newSVsv(hv_iterkeysv(ce)));
			}
			if (num_cols > 1)
				sortsv(AvARRAY(headers_av), num_cols, Perl_sv_cmp);
		}
		size_t num_headers = (size_t)(av_len(headers_av) + 1);
		const char **restrict header_row = safemalloc((num_headers + 1) * sizeof(char*));
		size_t h_idx = 0;
		if (inc_rownames) header_row[h_idx++] = "";
		for (size_t i = 0; i < num_headers; i++) {
			SV **restrict h_ptr = av_fetch(headers_av, (SSize_t)i, 0);
			header_row[h_idx++] = (h_ptr && SvOK(*h_ptr)) ? SvPV_nolen(*h_ptr) : "";
		}
		print_string_row(aTHX_ fh, header_row, h_idx, sep);
		safefree(header_row);
		const char **restrict row_data = safemalloc((num_headers + 1) * sizeof(char*));
		size_t d_idx = 0;
		// Give the single row a default numeric identifier if row names are on
		if (inc_rownames) row_data[d_idx++] = "1";
		for (size_t j = 0; j < num_headers; j++) {
			SV **restrict h_ptr = av_fetch(headers_av, (SSize_t)j, 0);
			SV *restrict h_sv = (h_ptr && SvOK(*h_ptr)) ? *h_ptr : NULL;
			// FIX (UTF-8/NUL safety): fetch by SV, not by raw bytes
			HE *restrict val_he = h_sv ? hv_fetch_ent(data_hv, h_sv, 0, 0) : NULL;
			SV *restrict val_sv = val_he ? HeVAL(val_he) : NULL;
			// FIX: a flat-hash cell holding a reference was stringified
			// (e.g. ARRAY(0x...)) instead of croaking like every other shape.
			if (val_sv && SvOK(val_sv)) {
				if (SvROK(val_sv)) {
					PerlIO_close(fh);
					safefree(row_data);
					if (headers_av) SvREFCNT_dec(headers_av);
					croak("write_table: Cannot write nested reference types to table\n");
				}
				row_data[d_idx++] = SvPV_nolen(val_sv);
			} else {
				row_data[d_idx++] = undef_val;
			}
		}
		print_string_row(aTHX_ fh, row_data, d_idx, sep);
		safefree(row_data);
	// ----- Hash of Arrays -----
	} else if (is_hoa) {
		HV *restrict data_hv = (HV*)data_ref;
		size_t max_rows = 0;
		hv_iterinit(data_hv);
		HE *restrict entry;
		while ((entry = hv_iternext(data_hv))) {
			AV *restrict arr = (AV*)SvRV(hv_iterval(data_hv, entry));
			size_t len = (size_t)(av_len(arr) + 1);
			if (len > max_rows) max_rows = len;
		}
		if (col_names_sv && SvOK(col_names_sv)) {
			AV *restrict c_av = (AV*)SvRV(col_names_sv);
			// FIX: size_t vs av_len() == -1 (empty col.names looped forever)
			for (SSize_t i = 0; i <= av_len(c_av); i++) {
				SV **restrict c = av_fetch(c_av, i, 0);
				if (c && SvOK(*c)) av_push(headers_av, newSVsv(*c));
			}
		} else {
			// FIX (UTF-8 safety): keep the key SVs (flags intact) and sort
			// them with sv_cmp instead of round-tripping through char*.
			unsigned int num_cols = hv_iterinit(data_hv);
			for (unsigned int i = 0; i < num_cols; i++) {
				HE *restrict ce = hv_iternext(data_hv);
				av_push(headers_av, newSVsv(hv_iterkeysv(ce)));
			}
			if (num_cols > 1)
				sortsv(AvARRAY(headers_av), num_cols, Perl_sv_cmp);
		}
		if (av_len(headers_av) < 0) {
			// FIX: this croak leaked the open filehandle and headers_av.
			PerlIO_close(fh);
			SvREFCNT_dec(headers_av);
			croak("Could not get headers in write_table");
		}
		if (inc_rownames && contains_nondigit(aTHX_ row_names_sv)) {
			rownames_col = SvPV_nolen(row_names_sv);
			AV *restrict filtered_headers = newAV();
			// FIX: size_t vs av_len() (same wrap as above if headers empty)
			for (SSize_t i = 0; i <= av_len(headers_av); i++) {
				SV **restrict h_ptr = av_fetch(headers_av, i, 0);
				if (!h_ptr || !*h_ptr) continue;
				SV *restrict h_sv = *h_ptr;
				// FIX (UTF-8 safety): sv_eq, not strcmp on raw bytes
				if (!sv_eq(h_sv, row_names_sv)) {
					av_push(filtered_headers, newSVsv(h_sv));
				}
			}
			SvREFCNT_dec(headers_av);
			headers_av = filtered_headers;
		}
		size_t num_headers = (size_t)(av_len(headers_av) + 1);
		const char **restrict header_row = safemalloc((num_headers + 1) * sizeof(char*));
		size_t h_idx = 0;
		if (inc_rownames) header_row[h_idx++] = "";
		for (size_t i = 0; i < num_headers; i++) {
			SV **restrict h_ptr = av_fetch(headers_av, (SSize_t)i, 0);
			header_row[h_idx++] = (h_ptr && SvOK(*h_ptr)) ? SvPV_nolen(*h_ptr) : "";
		}
		print_string_row(aTHX_ fh, header_row, h_idx, sep);
		safefree(header_row);
		const char **restrict row_data = safemalloc((num_headers + 1) * sizeof(char*));
		// FIX: numeric row labels used savepv() + safefree() every row; a
		// stack buffer reused per row does the same job with no allocation
		// (and removes the const-away cast in the old safefree call).
		char rn_buf[32];
		for (size_t i = 0; i < max_rows; i++) {
			size_t d_idx = 0;
			if (inc_rownames) {
				if (rownames_col) {
					// FIX (UTF-8 safety): fetch the row-name column by SV
					HE *restrict rn_arr_he = hv_fetch_ent(data_hv, row_names_sv, 0, 0);
					SV *restrict rn_arr_sv = rn_arr_he ? HeVAL(rn_arr_he) : NULL;
					if (rn_arr_sv && SvROK(rn_arr_sv)) {
						AV *restrict rn_arr = (AV*)SvRV(rn_arr_sv);
						SV **restrict rn_val_ptr = av_fetch(rn_arr, (SSize_t)i, 0);
						if (rn_val_ptr && SvOK(*rn_val_ptr)) {
							if (SvROK(*rn_val_ptr)) {
								PerlIO_close(fh);
								safefree(row_data);
								if (headers_av) SvREFCNT_dec(headers_av);
								croak("write_table: Cannot write nested reference types to table\n");
							}
							row_data[d_idx++] = SvPV_nolen(*rn_val_ptr);
						} else {
							row_data[d_idx++] = undef_val;
						}
					} else {
						row_data[d_idx++] = undef_val;
					}
				} else {
					snprintf(rn_buf, sizeof(rn_buf), "%lu", (unsigned long)(i + 1));
					row_data[d_idx++] = rn_buf;
				}
			}
			for (size_t j = 0; j < num_headers; j++) {
				SV **restrict h_ptr = av_fetch(headers_av, (SSize_t)j, 0);
				SV *restrict h_sv = (h_ptr && SvOK(*h_ptr)) ? *h_ptr : NULL;
				// FIX (UTF-8/NUL safety): fetch by SV, not by raw bytes
				HE *restrict arr_he = h_sv ? hv_fetch_ent(data_hv, h_sv, 0, 0) : NULL;
				SV *restrict arr_sv = arr_he ? HeVAL(arr_he) : NULL;
				if (arr_sv && SvROK(arr_sv)) {
					AV *restrict arr = (AV*)SvRV(arr_sv);
					SV **restrict cell_ptr = av_fetch(arr, (SSize_t)i, 0);
					if (cell_ptr && SvOK(*cell_ptr)) {
						if (SvROK(*cell_ptr)) {
							PerlIO_close(fh);
							safefree(row_data);
							if (headers_av) SvREFCNT_dec(headers_av);
							croak("write_table: Cannot write nested reference types to table\n");
						}
						row_data[d_idx++] = SvPV_nolen(*cell_ptr);
					} else {
						row_data[d_idx++] = undef_val;
					}
				} else {
					row_data[d_idx++] = undef_val;
				}
			}
			print_string_row(aTHX_ fh, row_data, d_idx, sep);
		}
		safefree(row_data);
	} else if (is_aoh) { // ----- Array of Hashes
		AV *restrict data_av = (AV*)data_ref;
		size_t num_rows = (size_t)(av_len(data_av) + 1);
		if (col_names_sv && SvOK(col_names_sv)) {
			AV *restrict c_av = (AV*)SvRV(col_names_sv);
			// FIX: size_t vs av_len() == -1 (empty col.names looped forever)
			for (SSize_t i = 0; i <= av_len(c_av); i++) {
				SV **restrict c = av_fetch(c_av, i, 0);
				if (c && SvOK(*c)) av_push(headers_av, newSVsv(*c));
			}
		} else {
			HV *restrict col_map = newHV();
			for (size_t i = 0; i < num_rows; i++) {
				SV **restrict row_ptr = av_fetch(data_av, (SSize_t)i, 0);
				if (row_ptr && SvROK(*row_ptr)) {
					HV *restrict row_hv = (HV*)SvRV(*row_ptr);
					hv_iterinit(row_hv);
					HE *restrict entry;
					while ((entry = hv_iternext(row_hv))) {
						hv_store_ent(col_map, hv_iterkeysv(entry), newSViv(1), 0);
					}
				}
			}
			unsigned num_cols = hv_iterinit(col_map);
			// FIX (UTF-8 safety): keep the key SVs (flags intact) and sort
			// them with sv_cmp instead of round-tripping through char*.
			for (unsigned int i = 0; i < num_cols; i++) {
				HE *restrict ce = hv_iternext(col_map);
				av_push(headers_av, newSVsv(hv_iterkeysv(ce)));
			}
			if (num_cols > 1)
				sortsv(AvARRAY(headers_av), num_cols, Perl_sv_cmp);
			SvREFCNT_dec(col_map);
		}
		if (inc_rownames && contains_nondigit(aTHX_ row_names_sv)) {
			rownames_col = SvPV_nolen(row_names_sv);
			AV *restrict filtered_headers = newAV();
			// FIX: size_t vs av_len() (same wrap as above if headers empty)
			for (SSize_t i = 0; i <= av_len(headers_av); i++) {
				SV **restrict h_ptr = av_fetch(headers_av, i, 0);
				if (!h_ptr || !*h_ptr) continue;
				SV *restrict h_sv = *h_ptr;
				// FIX (UTF-8 safety): sv_eq, not strcmp on raw bytes
				if (!sv_eq(h_sv, row_names_sv)) {
					av_push(filtered_headers, newSVsv(h_sv));
				}
			}
			SvREFCNT_dec(headers_av);
			headers_av = filtered_headers;
		}
		size_t num_headers = (size_t)(av_len(headers_av) + 1);
		const char **restrict header_row = safemalloc((num_headers + 1) * sizeof(char*));
		size_t h_idx = 0;
		if (inc_rownames) header_row[h_idx++] = "";
		for (size_t i = 0; i < num_headers; i++) {
			SV **restrict h_ptr = av_fetch(headers_av, (SSize_t)i, 0);
			header_row[h_idx++] = (h_ptr && SvOK(*h_ptr)) ? SvPV_nolen(*h_ptr) : "";
		}
		print_string_row(aTHX_ fh, header_row, h_idx, sep);
		safefree(header_row);
		const char **restrict row_data = safemalloc((num_headers + 1) * sizeof(char*));
		char rn_buf[32]; // FIX: replaces per-row savepv/safefree (see HoA)
		for (size_t i = 0; i < num_rows; i++) {
			size_t d_idx = 0;
			SV **restrict row_ptr = av_fetch(data_av, (SSize_t)i, 0);
			HV *restrict row_hv = (row_ptr && SvROK(*row_ptr)) ? (HV*)SvRV(*row_ptr) : NULL;
			if (inc_rownames) {
				if (rownames_col) {
					// FIX (UTF-8 safety): fetch the row-name cell by SV
					HE *restrict rn_he = row_hv ? hv_fetch_ent(row_hv, row_names_sv, 0, 0) : NULL;
					SV *restrict rn_sv = rn_he ? HeVAL(rn_he) : NULL;
					if (rn_sv && SvOK(rn_sv)) {
						if (SvROK(rn_sv)) {
							PerlIO_close(fh);
							safefree(row_data);
							if (headers_av) SvREFCNT_dec(headers_av);
							croak("write_table: Cannot write nested reference types to table\n");
						}
						row_data[d_idx++] = SvPV_nolen(rn_sv);
					} else {
						row_data[d_idx++] = undef_val;
					}
				} else {
					snprintf(rn_buf, sizeof(rn_buf), "%lu", (unsigned long)(i + 1));
					row_data[d_idx++] = rn_buf;
				}
			}
			for (size_t j = 0; j < num_headers; j++) {
				SV **restrict h_ptr = av_fetch(headers_av, (SSize_t)j, 0);
				SV *restrict h_sv = (h_ptr && SvOK(*h_ptr)) ? *h_ptr : NULL;
				// FIX (UTF-8/NUL safety): fetch by SV, not by raw bytes
				HE *restrict cell_he = (row_hv && h_sv) ? hv_fetch_ent(row_hv, h_sv, 0, 0) : NULL;
				SV *restrict cell_sv = cell_he ? HeVAL(cell_he) : NULL;
				if (cell_sv && SvOK(cell_sv)) {
					if (SvROK(cell_sv)) {
						PerlIO_close(fh);
						safefree(row_data);
						if (headers_av) SvREFCNT_dec(headers_av);
						croak("write_table: Cannot write nested reference types to table\n");
					}
					row_data[d_idx++] = SvPV_nolen(cell_sv);
				} else {
					row_data[d_idx++] = undef_val;
				}
			}
			print_string_row(aTHX_ fh, row_data, d_idx, sep);
		}
		safefree(row_data);
	}
	if (headers_av) SvREFCNT_dec(headers_av);
	if (rows_av) SvREFCNT_dec(rows_av);
	PerlIO_close(fh);
	XSRETURN_EMPTY;
}
