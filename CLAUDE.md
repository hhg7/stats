# Stats::LikeR — instructions for Claude

## Never edit `Changes`

Do not create, edit, append to, or revert `Changes` under any circumstances.
The release notes there are written by hand, by the maintainer, and are not
Claude's to touch — not even to add an entry for work Claude just did, and not
even when asked to "update the changelog" as part of a larger task.

This holds regardless of the tool: no `Edit`, no `Write`, no `sed -i`/`perl -pi`,
no `git checkout`/`git revert` that touches it, no patch that includes it.

When a change would normally warrant a release note, say so in the reply and
leave the wording to the maintainer. Do not draft it into the file.

## C types must match the value's real domain

Every variable and parameter in `LikeR.xs` gets the type that describes what
the value can actually be, not the type that happens to be convenient. Plain
`int` is not a default — reach for it only when the domain genuinely is
"signed, machine-word-ish, and not one of the cases below".

- **Two states only (0 or 1, yes/no, present/absent)** → `bool`, assigned
  `TRUE`/`FALSE` (perl's own, already in scope via `perl.h`; do not add
  `<stdbool.h>`). Never `int flag = 0`.
- **Cannot be negative** → an unsigned type, and give it the *widest*
  practical range so it cannot overflow on real input:
  - sizes, lengths, counts, and array indices → `size_t`
  - a length crossing the perl string API → `STRLEN`
  - an integer handed to or taken from perl → `UV` (and `IV` when it must be
    signed)
  - a fixed width that the algorithm depends on (hashes, bitsets,
    accumulators that must not wrap at 32 bits) → `uint64_t` / `uint32_t`
  - a byte → `unsigned char`
  Do not "save" a `size_t` by using `unsigned int` when the value tracks the
  size of something the caller controls.
- **A small, enumerated set of possibilities** (a mode, a kind, a strategy —
  say up to a handful of values) → `short int`, or `unsigned short int` when
  no sentinel is needed. Use signed `short int` when `-1` means "not yet
  decided", which is the existing idiom (`short int exact = -1;`,
  `short int na_mode = 0; // 0 = pairwise, 1 = omit, 2 = keep`). Always
  comment what each value means at the declaration.
- **A loop counter bounded by a small literal** → `unsigned short int`, as in
  `for (unsigned short int i = 0; i < 200; i++)`. A counter bounded by a
  runtime count is `size_t` instead.
- **Floating point** → `NV` throughout, never bare `double`, so the module
  keeps working on long-double and quadmath perls.

Widening or narrowing must never silently change a perl-visible signature or
conversion (`SvUV` vs `SvIV` vs `SvNV`, `%zu` vs `%" UVuf "`); fix the format
strings and the `Sv*` calls in the same edit. Also mind the mixed-sign
comparison warnings this can introduce — the build is expected to stay clean.

## Use `restrict` by default

Pointer parameters and pointer locals get `restrict` unless there is a real
reason they cannot. Apply it even when the compiler emits byte-identical
code: it documents the non-aliasing contract for the next reader, and it is
checked by review, not by the optimizer. `Makefile.PL` compiles with
`-std=gnu99`, so write bare `restrict` (not `__restrict__`).

Pair it with `const` when the pointee is not written, which is the existing
style: `const char *restrict alt`, `const NV *restrict lgR`,
`bool *restrict full`.

Leave `restrict` off — and say why in a short comment — when:

- two parameters may legitimately point at the same object or into the same
  buffer (in-place transforms, `src == dst` fast paths, overlapping
  slices/rows of one allocation);
- pointers are swapped or re-pointed at each other's targets during the
  function's lifetime;
- the pointee can be reached by another route inside the same call — through
  a struct field also passed separately, a recursive call, or a callback;
- it points at perl-managed memory that the perl API may alias behind your
  back: `SV`/`AV`/`HV` internals, a `PV` that may be COW or a shared hash
  key, or a buffer that can be reallocated by an intervening perl call;
- the pointer escapes to code you cannot audit.

Do not churn existing declarations solely to add or remove `restrict`; apply
this when writing new code or when already editing the function.

## Tests must come from R's and Python's own test suites

If a function has an equivalent in R or in Python (SciPy, NumPy, pandas,
statsmodels), tests for the `Stats::LikeR` version must be **taken from those
projects' own test suites and documented examples**, not invented here.
Invented cases only confirm that the code does what it currently does; the
references' cases are what pin it to the behavior it is supposed to copy.

Where to look — both suites are available locally, so read them rather than
recalling them:

- R 4.6.1 source: `/home/con/Scripts/r-source` —
  `tests/reg-tests-1{a,b,c,d,e}.R` and `tests/reg-tests-2.R` for regression
  cases, `src/library/stats/man/*.Rd` for documented examples, and
  `tests/Examples/stats-Ex.Rout.save` for those examples' pinned output.
- SciPy 1.18.0: `/home/con/.pyenv/versions/3.14.2/lib/python3.14/site-packages/scipy/stats/tests/`
  — most `Test<Name>` classes are in `test_stats.py` and the hypothesis-test
  classes in `test_morestats.py`, but the split moves between releases, so grep
  the whole directory rather than those two files: 1.18.0 pulled
  `TestSpearmanRho`, `TestTheilslopes` and `TestSiegelslopes` out into
  `test_correlation.py` while leaving `TestPearsonr`, `TestCorrSpearmanr` and
  `TestKendallTau` behind in `test_stats.py`, and `test_rank.py`,
  `test_quantile.py` and `test_hypotests.py` are each worth checking too.
  `data/` holds R-generated reference corpora that are directly reusable as
  Perl test data.
- NumPy 2.5.2: same `site-packages`, tests under `numpy/_core/tests/` and
  `numpy/lib/tests/` (each with its own `data/`).

Check **both** before deciding a function is well covered: coverage is
lopsided, and the thinner suite is not always the one you would guess (R has
a single `binom.test` case in all of `reg-tests-*.R`, while SciPy pins dozens
for the same function, themselves annotated "results from R").

### Freeze the expected values; never depend on R or Python at run time

A generated test must pass on a machine with no R, no `python3`, no NumPy and
no SciPy. That means:

- Expected values are **hardcoded literals in the `.t` file** (or in a plain
  Perl data table inside it), at full precision — `options(digits=17)` on the
  R side, `repr()` on the Python side.
- No `Rscript`, no `python3`, no `qx{}`, no `system()`, and no
  `plan skip_all` / `SKIP` block conditioned on a reference implementation
  being installed. A cross-validation test that skips is a test that silently
  never runs on a CPAN smoker.
- When a generator script is worth keeping, commit it next to the test and
  say in a comment how to re-run it, as `t/bedroc.python.py` does for
  `t/bedroc.python.t`. The script regenerates the frozen table; the test
  still never calls it.
- Record the provenance of every expected value in the file's header comment:
  which reference, which version, which file, which test function or man
  page, and any upstream annotation about where the number itself came from.
  `t/binom_test.R.scipy.t` is the model to follow, including its naming
  (`t/<fn>.R.scipy.t` for the cross-validation file, kept separate from the
  hand-written `t/<fn>.t`).

### Maximum coverage

Aim to exercise every path, not one representative call:

- Cross the whole argument space the references cross — every alternative /
  tail / method / correction switch, and every combination of them, rather
  than the default.
- Cover the edges the references cover: n = 0 and 1, empty and single-element
  input, all-identical values, ties, zero variance, zero cells, exact
  boundaries, values that underflow or overflow, `NaN`/`Inf`, and every
  `NA`-handling mode.
- Cover the Perl-side surface too: each accepted call form (hashref, arrayref,
  named arguments), each documented return field, every `croak` message and
  argument-validation path, and leak checks where the sibling tests have them.
- Choose a tolerance per quantity, and justify it in a comment against the
  worst observed disagreement, leaving headroom for long-double and quadmath
  builds. Do **not** widen a tolerance to make a failure go away.
- Gate genuinely slow cases behind `EXTENDED_TESTING`/`AUTHOR_TESTING`
  rather than dropping them.

When the XS and a reference disagree, do not silently pick a winner and do
not paper over it: get a third opinion from `mpmath` at `mp.dps = 60`,
solving the *defining* equation by bisection rather than calling a library
inverse, and record any remaining divergence explicitly in the test file so
that changing it later is a deliberate act.

## Every change must hold across the whole support matrix

A modification is not finished when it works on the default perl. It has to
hold for every installed perl, all three NV widths, Windows, and perl 5.10.

### All installed perls

`./test.all.perls.pl` builds and tests against every perlbrew perl; run it (or
at minimum the affected versions with `-p`) for anything touching `LikeR.xs`.
The local matrix under `/home/con/perl5/perlbrew/perls/`, by `$Config{nvtype}`:

- `perl-5.44.0` (default), `perl-5.42.3`, `perl-5.10.1` — `double`, NVgf `"g"`
- `perl-5.12.5` — `long double`, NVgf `"Lg"` (archname `x86_64-linux-ld`)
- `5.44.0-quadmath` — `__float128`, NVgf `"Qg"`

All three NV widths can be compile-checked without reconfiguring the tree
(which would clobber the current `Makefile`) by generating the `.c` and
compiling it against each perl's `CORE` directly. Do that as the fast check;
it is not a substitute for actually running the suite.

### Long double and quadmath

- Floating point is `NV`, never bare `double`; use `NV_INF`, `NV_NAN`,
  `NV_EPSILON`, `NV_MAX` for the constants.
- **Never call libm bare.** `sqrt(x)`, `log(x)`, `lgamma(x)` and friends take a
  `double`, so on a long-double or `__float128` build they silently round the
  `NV` down to 53 bits, compute there, and convert back — no warning, no build
  failure, just a less accurate answer than the perl running it. Every
  NV-valued call goes through the `nv_*` macros defined below the includes in
  `LikeR.xs` (`nv_sqrt`, `nv_log`, `nv_lgamma`, …), which paste on the `l` or
  `q` suffix for the build's width. Adding a new libm function means adding it
  to that block *and* to the link probe in `Makefile.PL`, not calling it
  directly.
- Float classification goes through `nv_isnan`/`nv_isinf`/`nv_isfinite`,
  defined next to the `nv_*` libm block, and through nothing else. Not the bare
  C99 macros: where a platform does not provide the type-generic versions,
  `isfinite()` is a plain `double` function and a large-but-finite `NV`
  narrowed into it is reported as infinite. And **not** perl's
  `Perl_isnan`/`Perl_isinf`/`Perl_isfinite`, which is what `LikeR.xs` used up
  to 0.298: on every perl before 5.22 they can route through a
  `Perl_fp_class()` block in `perl.h` that has never compiled — the macro is
  spelled with an empty parameter list and compares against `FP_CLASS_*` names
  no `<ieeefp.h>` defines. Configure only reaches that block where it does not
  find `isinf()`, so it is dead code on Linux and glibc and live on
  illumos/Solaris, which is exactly how it reached CPAN in 0.298. Those perls
  also define `Perl_isinf(x)` as `((x)==NV_INF)`, missing `-Inf`.
- Format NVs with perl's `my_snprintf`/`my_sprintf` and `NVgf`/`NVff`, never
  the plain C-library `snprintf` with a hand-written `%g`. `%Qg` happens to
  work with glibc's `snprintf` only because libquadmath registers the `Q`
  specifier at load time — a glibc extension that does not exist elsewhere.
- Do not assume a literal or a constant is exactly representable, and do not
  hardcode a `double`-sized epsilon; scale tolerances off `NV_EPSILON`.
- Test tolerances need headroom for the wider widths, and any case whose cost
  scales with NV width (long-double or `__float128` libm in an inner loop) is
  a candidate for an `EXTENDED_TESTING` gate rather than a slow smoker.

### Windows, Solaris and every BSD

There is no local Windows, Solaris, illumos, FreeBSD, OpenBSD, NetBSD or
DragonFly perl, so all of this is discipline applied while writing, not
something a test run here will catch. The module is expected to build and pass
on all of them.

Common to all of these targets:

- Allocate with `Newx`/`Newxz`/`Renew`/`Safefree`, not `malloc`/`free`, so the
  module uses whatever allocator its perl was built with.
- Reach for the perl API before libc whenever perl has an equivalent
  (`my_snprintf`, `sortsv`, `strEQ`/`foldEQ`, `PerlIO`, `PerlProc_*`): it is
  the layer that has already been ported to every one of these platforms.
- Nothing may depend on a GNU extension: no VLAs, no nested functions, no
  statement expressions, no `typeof`, no `__builtin_*` without a fallback, no
  zero-length arrays, no case ranges, no `void *` arithmetic. `LikeR.xs` must
  keep compiling under strict `-std=c99`, which is the closest local proxy for
  a vendor compiler — check with `-std=c99` (not `gnu99`) when in doubt.
- C99 `restrict` is spelled bare in this file, which not every target accepts
  (MSVC wants `__restrict` outside `/std:c11`, older gcc `__restrict__`, a C89
  compiler has nothing). The `#if !defined(__cplusplus) && !defined(restrict)`
  block below the includes maps or defines it away; keep that block, and if you
  add `restrict` to a translation unit that lacks one, add the guard too.
- The C99 flag itself is per-vendor and is probed, not guessed, by
  `Makefile.PL` — `-std=gnu99`/`-std=c99` for gcc and clang, `-xc99=all` for
  Oracle Studio, `-qlanglvl=extc99` for AIX `xlc`, `-AC99` for HP-UX, nothing
  for MSVC. `$Config{cc}` is plain `cc` on Solaris and HP-UX, so it can never
  be pattern-matched for a gcc-only flag. That logic lives in `dist.ini`'s
  `[MakeMaker::Awesome]` header (the root `Makefile.PL` is generated from it and
  mirrors it); change it in `dist.ini`, not only in `Makefile.PL`.

Windows specifically:

- No POSIX-only headers or functions: no `<unistd.h>`, `<sys/time.h>`, or
  `<strings.h>`; no `strcasecmp`/`strncasecmp` (there is `str_ieq_ascii()` in
  `LikeR.xs` for the case-insensitive keyword compares), no `fork`, no
  hardcoded `/` in paths.
- Printf lengths must be perl's: `%" UVuf "`, `%" IVdf "`, `%" NVgf "`,
  `%" UVxf "`. `%zu`, `%zd`, `%lld` and `%llu` are not portable to MSVC's
  older CRT — cast to `UV`/`IV` and use the perl macro.

Solaris, illumos and the BSDs specifically:

- Keep the `_GNU_SOURCE` / `__EXTENSIONS__` block at the top of `LikeR.xs`
  intact and add to it rather than around it: `__EXTENSIONS__` is what exposes
  the POSIX and XPG symbols on Solaris once `-std=c99`-style strictness is on,
  and dropping it silently loses declarations.
- Do not reach for glibc-only or BSD-only libc: no `qsort_r` (its argument
  order differs between glibc, the BSDs and Solaris — sort with perl's
  `sortsv`/`sortsv_flags` or a plain `qsort` with a file-static comparator),
  no `memmem`, `strchrnul`, `strndup`, `asprintf`, `getline`, `fmemopen`,
  `reallocarray`, `explicit_bzero`, `arc4random`, `exp10` or `sincos`.
- Stay inside C99 `<math.h>` for float math (`lgamma`, `log2`, `round`,
  `trunc`, `fmax` are all there); anything beyond it — `lgamma_r`, `j0`,
  `tgammal`, `lgammaq` — is not uniformly available and needs a guard plus a
  fallback. OpenBSD in particular is thin on the long-double variants.
- `/dev/urandom` exists on all of these, but a randomness path must still keep
  its documented fallback for the ones where opening it can fail.
- Assume nothing about a struct's field order, `char`'s signedness, or
  alignment: SPARC and some ARM/BSD combinations are strict about unaligned
  access, so never cast a `char *` to a wider pointer type and dereference it.

### Back-compatible to perl 5.10

`dist.ini` declares `perl = 5.010` and `lib/Stats/LikeR.pm` says
`require 5.010`; both the Perl and the C sides must honor that.

- Perl code (module and tests) is limited to 5.10 syntax: no signatures, no
  postfix dereference, no `say` without importing the feature, no `isa`, no
  chained comparisons, no `//=`-era-only assumptions beyond what 5.10 has.
  Tests follow the existing header — `require 5.010; use strict; use warnings
  FATAL => 'all';` — not `use 5.044`. The author-only helper scripts in the
  repo root (`test.all.perls.pl`, `xs.check.pl`) are exempt; they are not
  shipped and may use modern perl.
- Any perl C API newer than 5.10 must come through `ppport.h`, and
  `ppport.h` must be regenerated when a new one is used. If `ppport.h` does
  not cover it, write the fallback rather than raising the minimum.
- Do not use API that was merely *renamed* later, and do not rely on
  behavior that changed after 5.10 (hash iteration order, `SvPV` return
  guarantees, COW string semantics) without checking on `perl-5.10.1`.
- Raising `MIN_PERL_VERSION` is a maintainer decision, not a way to make a
  build error go away. If a change truly needs a newer perl, say so in the
  reply and stop.
