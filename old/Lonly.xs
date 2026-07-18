/*
 * Retired XSUB: the original two-argument Lonly().
 *
 * Removed from LikeR.xs when the `Lonly` name was reassigned to the former
 * `get_unique` (an @-prototype, N-ref "values only in the first list" filter
 * backed by set_multiplicity()). For two array refs the new Lonly() returns the
 * same result as this old one, so callers passing exactly two lists are
 * unaffected.
 *
 * This XSUB is backed by the shared set_difference() engine, which REMAINS live
 * in LikeR.xs because Ronly() still uses it. A copy of that engine is included
 * below for reference so this archive reads standalone; it is not compiled from
 * here.
 */

void Lonly(...)
	PROTOTYPE: $$
	PPCODE:
		if (items != 2)
			croak("Lonly needs exactly 2 array refs (got %" UVuf ")", (UV)items);
		if (!(SvROK(ST(0)) && SvTYPE(SvRV(ST(0))) == SVt_PVAV))
			croak("Lonly: first (left) argument is not an array reference");
		if (!(SvROK(ST(1)) && SvTYPE(SvRV(ST(1))) == SVt_PVAV))
			croak("Lonly: second (right) argument is not an array reference");
		/* keep = left (ST0), other = right (ST1) */
		SP = set_difference(aTHX_ SP, ST(0), ST(1), "left", "right", "Lonly", GIMME_V);

/* ------------------------------------------------------------------------- *
 * Supporting engine (reference copy; the live definition is in LikeR.xs and
 * is still used by Ronly()):
 *
 * Backs Lonly() and Ronly(): the values in `keep` (deduped, first-seen order)
 * that do not occur in `other`. Ronly is just Lonly with the two arrays
 * swapped, so both XSUBs call this.
 * ------------------------------------------------------------------------- */
static SV** set_difference(pTHX_ SV **sp, SV *keep_sv, SV *other_sv,
                           const char *keep_side, const char *other_side,
                           const char *name, int gimme) {
	HV *restrict inOther = (HV*)sv_2mortal((SV*)newHV());
	HV *restrict seen    = (HV*)sv_2mortal((SV*)newHV());
	AV *restrict other   = (AV*)SvRV(other_sv);
	AV *restrict keep    = (AV*)SvRV(keep_sv);
	size_t olen = (size_t)(av_len(other) + 1);
	size_t klen_a, n = 0;
	/* keep_side/other_side name the physical (left/right) position of each
	 * array so the croak matches how the value was passed, regardless of the
	 * keep/subtract role (Ronly swaps the roles but not the positions). */
	for (size_t j = 0; j < olen; j++) {
		SV **restrict tv = av_fetch(other, j, 0);
		STRLEN klen; const char *restrict key; I32 hklen;
		if (!(tv && SvOK(*tv)))
			croak("%s: undefined value in %s array ref at index %" UVuf, name, other_side, (UV)j);
		key   = SvPV(*tv, klen);
		hklen = SvUTF8(*tv) ? -(I32)klen : (I32)klen;
		(void)hv_store(inOther, key, hklen, &PL_sv_undef, 0);
	}
	klen_a = (size_t)(av_len(keep) + 1);
	for (size_t j = 0; j < klen_a; j++) {
		SV **restrict tv = av_fetch(keep, j, 0);
		STRLEN klen; const char *restrict key; I32 hklen;
		if (!(tv && SvOK(*tv)))
			croak("%s: undefined value in %s array ref at index %" UVuf, name, keep_side, (UV)j);
		key   = SvPV(*tv, klen);
		hklen = SvUTF8(*tv) ? -(I32)klen : (I32)klen;
		if (hv_exists(inOther, key, hklen)) continue;   /* present in other set   */
		if (hv_exists(seen,    key, hklen)) continue;   /* dedup within keep set  */
		(void)hv_store(seen, key, hklen, &PL_sv_undef, 0);
		n++;
		if (gimme != G_SCALAR) XPUSHs(sv_2mortal(newSVsv(*tv)));
	}
	if (gimme == G_SCALAR) XPUSHs(sv_2mortal(newSVuv(n)));
	return sp;
}
