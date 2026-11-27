# Stats::LikeR — instructions for Claude

## Editing `Changes`

Claude may edit `Changes`; the maintainer lifted the old prohibition on
2026-09-04. It is now the only copy of the release notes, nothing regenerates
it, and CPAN and PAUSE parse it, so a bad edit is a lost release note.

- Add a release by prepending a section to `Changes`, never to README.md
  (which no longer carries the notes). Never re-word a release that has
  shipped.
- A version line is `<version> <date> <tz>`, as in `0.315 2026-09-04 CDT`. A
  missing date or an unparsable version fails `changes_file_ok()` in
  `md2pod.pl`, which only checks the file now (up to 0.315 it regenerated it
  from README.md). After a hand-edit, load the file through `CPAN::Changes`
  and confirm it parses.
- Group entries under `[Bracketed headings]`, wrap at the width of the
  surrounding releases, and say what changed and why. The 0.315 and 0.314
  sections are the model.
- `$VERSION` in `lib/Stats/LikeR.pm` is the version *in progress*, not a
  released one; local tarballs and commits named after a version are not
  proof of release either. Check what has shipped first. Work done while
  `$VERSION` is unreleased goes under that version's section; do not open the
  next one.
- `release-notes.staged` is an untracked working copy of the whole file for
  staging notes, and can be ahead of `Changes`. Compare the two before
  assuming either is current; re-sync or delete it once its notes land.

The prohibition comes back when `.claude/hooks/no-edit-changes.sh` is
registered again as a `PreToolUse` hook in `.claude/settings.json` (emptied
2026-09-04); while registered it blocks every write to `Changes` whatever this
file says. It blocks `Edit`/`Write` on a basename `Changes` and any `Bash`
command that names `Changes` and might mutate it, including false positives
with a `>` before the name, such as
`perl -MCPAN::Changes -e 'CPAN::Changes->load("Changes")'`.

## Patch leftovers never get committed

After applying a patch by hand, delete its `.orig` and `.rej` files before
committing. `[PruneFiles]` ships everything under `lib/` and `t/`, and
`.gitignore` and `dzil.sh` cannot keep out a file that is already tracked:
`lib/Stats/LikeR.pm.orig` shipped in 0.3415 and 0.3416, and
`lib/Stats/LikeR.pm.rej` in 0.3422, where `make` installed it beside the
module.

## Comments in `LikeR.xs` use `//` when they fit on one line

A comment that fits on one line is a `//` comment, standalone or trailing,
with a space after the `//`:

    short int exact = -1; // -1 = decide from the sample sizes

A longer comment keeps the `/* ... */` block form: the text opens on the `/*`
line, continuation lines start at column 0 with no `*` gutter, and `*/` closes
the last line of text.

A line that ends in `\` (a macro continuation) keeps `/* ... */` even when the
comment is short, because a `//` there swallows the backslash and ends the
macro.

## C types must match the value's real domain

Every variable and parameter in `LikeR.xs` gets the type that describes what
the value can actually be. Plain `int` is only for a value that is genuinely
signed, machine-word-ish, and none of the cases below.

- **Two states** → `bool`, assigned `TRUE`/`FALSE` (perl's, from `perl.h`; do
  not add `<stdbool.h>`). Never `int flag = 0`.
- **Cannot be negative** → unsigned, at the widest practical range:
  - sizes, lengths, counts, array indices → `size_t`
  - a length crossing the perl string API → `STRLEN`
  - an integer to or from perl → `UV` (`IV` when it must be signed)
  - a fixed width the algorithm depends on (hashes, bitsets, accumulators
    that must not wrap at 32 bits) → `uint64_t` / `uint32_t`
  - a byte → `unsigned char`
  Never `unsigned int` for something whose size the caller controls.
- **A small enumerated set** (mode, kind, strategy) → `short int`, or
  `unsigned short int` when no sentinel is needed; signed when `-1` means
  "not yet decided" (`short int exact = -1;`). Comment every value at the
  declaration: `short int na_mode = 0; // 0 = pairwise, 1 = omit, 2 = keep`.
- **A loop counter bounded by a small literal** → `unsigned short int`; by a
  runtime count → `size_t`.
- **Floating point** → `NV`, never bare `double`.

A type change must not silently change a perl-visible signature or conversion:
fix the `Sv*` calls (`SvUV`/`SvIV`/`SvNV`) and format strings in the same
edit, and keep the build free of mixed-sign comparison warnings.

## Use `restrict` by default

Pointer parameters and locals get bare `restrict` (`Makefile.PL` compiles with
`-std=gnu99`), paired with `const` when the pointee is not written:
`const char *restrict alt`, `bool *restrict full`. Apply it even when the code
generated is identical; it documents the non-aliasing contract for review.

Leave it off, with a short comment saying why, when:

- two parameters may point at the same object or into one buffer (in-place
  transforms, `src == dst` paths, overlapping rows of one allocation);
- pointers are swapped or re-pointed at each other's targets;
- the pointee is reachable another way in the same call (a struct field also
  passed separately, recursion, a callback);
- it points at perl-managed memory the API may alias (`SV`/`AV`/`HV`
  internals, a COW `PV` or shared hash key, a buffer an intervening perl call
  may reallocate);
- the pointer escapes to code you cannot audit.

Do not churn existing declarations just to add or remove it; apply it to new
code and to functions already being edited.

## Loop variables live inside their loop

A loop counter, and anything a loop uses that nothing after it reads, is
declared in the loop (`for (size_t j = 0; ...)`, temporaries in the innermost
block that uses them), so a stale or end value cannot be read by mistake.

- Each loop declares its own counter; consecutive loops do not share one.
- The type rules still apply: `size_t` for a runtime bound,
  `unsigned short int` for a small literal one, `SSize_t` when the bound is
  `AvFILLp()`/`av_len()` and the comparison would otherwise mix signs.
- **Exception:** a value read after the loop (where a scan stopped, a found
  index, a count, as `if (k == first)` in `xlsx_xml_uncat()`) is declared
  before it, with a short comment saying it is read afterwards.
- `for` init declarations are C99, which `Makefile.PL` requires and MSVC has
  accepted since VS 2013.
- Perl code follows the same rule: `for my $i (...)`, `my` inside the body.

Do not churn functions just to move declarations; apply this to new code and
functions already being edited.

## Tests

Tests for any function with an R or Python equivalent come from those
projects' own test suites, frozen as literals. The full rules (sources,
freezing, coverage, tolerances) are in `.claude/rules/tests.md`, which loads
when you work on files under `t/`; read it before writing a test.

## Every change must hold across the whole support matrix

That means every installed perl, all three NV widths, both IV widths, Windows,
and perl 5.10.

### Installed perls

`./test.all.perls.pl` builds and tests against every perlbrew perl (`-p` for
a subset); run it for anything touching `LikeR.xs`. `--list` prints the
matrix with each perl's NV and IV width. Three perls are there for one
reason each: `perl-5.32.1` for that release's own bugs, `5.16.3-thr-ld` as
the only threaded perl older than 5.20, and `5.44.0-debug` for perl's
internal `assert()`s (see below).

As a fast check, compile-check all three NV widths without reconfiguring the
tree (which would clobber the `Makefile`) by generating the `.c` and compiling
it against each perl's `CORE`. It does not replace running the suite.

### A perl's own bugs show only on that perl

Varying widths and threading cannot catch a bug in one release's internals.
0.3214 failed 14 leak subtests on 5.32.1 only, because 5.32's regcomp leaks
the `\s`/`\w` inversion lists of every pattern compiled at run time, and
`_xml_attr()`/`_xml_start_tags()` interpolated a name into their patterns.

- In a function a leak test covers, never interpolate into a match that runs
  more than once per name. Compile it once into a `qr//` in a lexical cache,
  evaluated only on a cache miss (which also keeps the 5.10.0 `pp_qr()` leak
  below out of the warmed check).
- When a smoker fails on a series with no perl here, add one and reproduce
  there before changing anything; `./test.all.perls.pl --list` names the
  missing series and prints how to install one.

### A DEBUGGING perl asserts what a normal one allows

A perl configured with `-DDEBUGGING` compiles in its internal `assert()`s; a
normal perl silently accepts what they reject. Several smokers run only
DEBUGGING perls, so such a bug fails every one of their reports and none of
anyone else's. All 11 FAILs of 0.3422 were this: `S_av_push_own()` handed
`av_push()` the NULL `S_push_cell()` stores for a column `keep_cols` or
`drop_cols` leaves out, and `Perl_av_push` aborted with "Assertion `val'
failed". Those smokers' perls were also long double or quadmath, which was
the testers' configuration, not the cause: check ccflags for `-DDEBUGGING`
before trusting an NV-width pattern in the reports.

- Never pass NULL to `av_push`, `av_store`, `hv_store` or any other API that
  takes a value. A NULL slot is written directly into an array this code
  owns and has checked has no magic, after `av_extend()` if it is full, as
  `S_av_push_own()` does.
- `5.44.0-debug` is in the matrix for this. To check another API's assertion
  without it, generate `LikeR.c` in a scratch copy of the tree, redefine the
  macro after `#include "XSUB.h"` to abort on the asserted condition, build,
  and run the whole suite.

### The IV is not always 64 bits

With `use64bitint=undef` (every 32-bit smoker, armv6l included) `IV_MAX` is
2147483647, independently of the NV: a value exact as an `NV` can be far
outside `IV` range. So:

- Gate every `NV` → `IV`/`UV` cast on `IV_MAX`/`UV_MAX` as well as on 2\*\*53.
- Show that `IV` arithmetic on caller-supplied values cannot overflow at 32
  bits.

0.314 gated `seq()`'s integer path on 2\*\*53 alone, and `seq(1e15, 1e15+200,
2)` returned 2147483647 on armv6l. `t/seq.R.t` covers both branches, keying
expectations off `$Config{ivsize}`, and fails 30 subtests on `5.44.0-i686`
with the old gate.

Such a perl cannot be configured on x86-64 (Configure takes `ivsize` from
`longsize` and refuses `ivsize < ptrsize`); it must be a real 32-bit build,
whose recipe `./test.all.perls.pl` prints when the matrix lacks one. The
comment above that recipe in the script explains its two non-obvious flags
(`-Dplibpth`, `-Dlibswanted`).

### A size_t is not always as wide as an IV

On a 32-bit `-Duse64bitint` perl the `IV` is 8 bytes and `size_t` 4, so a
caller's integer can fit an `IV` and truncate at 2\*\*32 in a `size_t`.
`size_t` is right for counts of things in memory; a quantity that is not one,
such as `rbinom()`'s `size` (a trial count), is a `UV`/`IV`. 0.3414 capped it
with `(size_t)IV_MAX`, refused `size => 1e15`, and failed `t/rbinom.dist.t`
85–86 on `i686-linux-64int`; `5.44.0-i686-64int` reproduces it (built with the
iv4 recipe, `-Duse64bitint` instead of `-Uuse64bitint`).

### Whether an allocation fails is the machine's decision

- No test may depend on the machine refusing an allocation: an overcommitting
  system grants it and kills the process on first touch. Test an
  out-of-memory guard only with a request larger than the address space,
  such as `runif(2**48)` (2\*\*48 is `SV_COUNT_MAX`; the values need 2\*\*53.3
  bytes), and gate anything smaller behind `EXTENDED_TESTING`/`AUTHOR_TESTING`
  with a comment saying why.
- A `malloc()` only compared with NULL and freed may be deleted by the
  compiler, so the probe answers yes to any size. Write to the block through
  a volatile lvalue before freeing it, as `alloc_probe_ok()` does.

0.3419 failed 11 of 36 reports in `t/alloc.guards.t` for both reasons ("Out of
memory" after a probe on gcc 15.2/16.1, Hurd and Android clang; a granted
16 GB `bw_sj(nb => 2**31 - 1)` killed on NetBSD and DragonFly). The compiler
deletion fits the reports but was not observed locally. Forcing the probe to
`TRUE` in a 0.3419 build reproduces the Android report, and is how to check
that a new guard's test really exercises the probe.

### The C library decides how a tie rounds

perl hands `%f`, `%e` and `%g` to the C library, and the libraries round a
value exactly halfway between two results differently: glibc, UCRT,
libquadmath and R to the even digit, the MSVCRT under Strawberry Perl 5.24
away from zero. Only perl's own `%.0f` rounds the same everywhere.

- Never format a number a user sees at a fixed precision with `sprintf`
  (`'%.1f'`, `'%.*f'`, `'%.4g'`) or `my_snprintf`. Use `_format_nv($style,
  $prec, $x)` in perl, or `sv_set_nv_even()` in `LikeR.xs`, which round an
  exact tie to even themselves and leave every other value to the library.
  `%g` precisions above 13, and ties with more than 17 significant digits,
  are still the library's.
- A test that compares formatted numbers includes ties whose even neighbor
  is below them (0.125 at 2 places, 56.25 at 1), since glibc gets these
  right and so passes either way. Before trusting the test, run it under an
  `LD_PRELOAD` shim that rounds ties away from zero.

0.3431's `table_one()` printed 18 of 32 as "18 (56.3%)" on Strawberry 5.24
and "18 (56.2%)" elsewhere, failing `t/table_one.R.datasets.t` 120. The shim
reproduced it on Linux, failing that subtest and no other in the suite.
`summary()` and `write_table(tex_format => 1)` had the same bug, unreported.
`t/format.ties.R.t` is the regression test.

### A smoker can be 300 times slower

The armv6l-linux smoker is a Raspberry Pi. It runs the suite 48 times slower
than Strawberry, and an 8001-term `lm()` fit 300 times slower than x86-64.

- A wall-clock bound in a test covers only the code the test is about. Keep
  out of the timed region any work that grows faster than linearly, or needs
  hundreds of MB, by nature: a fit, an X'X, a sort of the input.
- State the bound's headroom against the time measured here, and check that
  it still holds at 300 times slower.

0.3431 failed `t/formula.parser.R.t` 119 there at 150 s against a 120 s
bound. The subtest checked that the formula parse is linear, but it timed
the whole `lm()`, whose X'X and its factor are 8001 x 8001 (760 MB peak,
0.5 s here). It now times the parse alone, through the "object 'v1' not
found" path, which takes 5 ms here.

### Threads change which libc function perl calls

Under `USE_REENTRANT_API`, `reentr.h` maps many libc interfaces (`drand48`,
`random`, `strerror`, `readdir`, `localtime`, `setlocale`, `getpw*`/`getgr*`,
…) onto their `_r` forms over `PL_reentrant_buffer` — but only under
`PERL_REENTR_API == 1`, i.e. for `PERL_CORE`/`PERL_EXT`, so perl's core gets
the `_r` form and `LikeR.xs` the plain one (`perlxs`, "Thread-aware system
interfaces", warns against mixing them). Never call a libc function on
`reentr.h`'s list from `LikeR.xs` without first checking what perl's core does
with it. 0.316's `Drand01()` (libc `drand48()` before 5.20) ignored
`srand($seed)` on threaded perls; the block below `AUTO_SEED_PRNG()` in
`LikeR.xs` is the fix, `t/srand.stream.t` the regression test, and
`5.16.3-thr-ld` the perl that fails without it (`perl-5.42.3` is too new to
show it).

### Long double and quadmath

- Floating point is `NV`, with `NV_INF`, `NV_NAN`, `NV_EPSILON`, `NV_MAX`.
- **Never call libm bare**: `sqrt()`, `log()`, `lgamma()` and the rest take a
  `double` and silently round a wider `NV`. Use the `nv_*` macros below the
  includes in `LikeR.xs`, which add the `l`/`q` suffix. A new libm function
  goes into that block *and* the link probe in `Makefile.PL`.
- **A libm function that links is not necessarily correct.** newlib's
  `fma()`, Cygwin's libm, is `x * y + z` with two roundings, so 0.3428's
  `ib_lambda()`, which needed one rounding, left `pbinom()` 2.5e-10 off at
  n = 1.8e13 there (`t/numeric.fixes.t` 38). Code whose accuracy rests on a
  libm guarantee beyond "close" (one rounding, exactness, correct rounding)
  checks it at load time with volatile operands, as `nv_fma_probe()` does in
  `BOOT`, and keeps a software fallback (`nv_two_prod()`). To test the
  fallback, redefine the function badly in a scratch build and run the suite.
- Classify floats only with `nv_isnan`/`nv_isinf`/`nv_isfinite`. Not the bare
  C99 macros (where not type-generic, `isfinite()` takes a `double` and calls
  a large finite `NV` infinite), and not `Perl_isnan`/`Perl_isinf`/
  `Perl_isfinite`: before 5.22 they can route through a `Perl_fp_class()`
  block in `perl.h` that has never compiled, live where `isinf()` is not
  found (illumos/Solaris, which is how 0.298 broke), and `Perl_isinf(x)` is
  `((x)==NV_INF)`, missing `-Inf`.
- Format NVs with `my_snprintf`/`my_sprintf` and `NVgf`/`NVff`, never libc
  `snprintf` with a hand-written `%g`; `%Qg` works in glibc only through a
  libquadmath extension. A number a user sees at a fixed precision goes
  through `sv_set_nv_even()` instead (see "The C library decides how a tie
  rounds").
- Assume no literal is exactly representable; scale tolerances off
  `NV_EPSILON`, never a hardcoded `double` epsilon.
- Give test tolerances headroom for wider NVs, and consider an
  `EXTENDED_TESTING` gate for cases whose cost scales with NV width.

### Windows, Solaris and every BSD

None of these is available locally, so this is discipline while writing; the
module must still build and pass on all of them.

- Allocate with `Newx`/`Newxz`/`Renew`/`Safefree`, not `malloc`/`free`.
- Prefer the perl API to libc wherever it has an equivalent (`my_snprintf`,
  `sortsv`, `strEQ`/`foldEQ`, `PerlIO`, `PerlProc_*`).
- No GNU extensions: VLAs, nested functions, statement expressions,
  `typeof`, `__builtin_*` without a fallback, zero-length arrays, case
  ranges, `void *` arithmetic. `LikeR.xs` must compile under strict
  `-std=c99` (not `gnu99`), the local proxy for a vendor compiler.
- Keep the `#if !defined(__cplusplus) && !defined(restrict)` block below the
  includes, which maps or removes `restrict` for compilers that lack it; add
  the same guard to any other translation unit that uses `restrict`.
- The C99 flag is probed per vendor, not guessed (`$Config{cc}` is plain `cc`
  on Solaris and HP-UX): `-std=gnu99`/`-std=c99`, `-xc99=all`,
  `-qlanglvl=extc99`, `-AC99`, nothing for MSVC. That logic lives in
  `dist.ini`'s `[MakeMaker::Awesome]` header, which generates `Makefile.PL`;
  change it there.

Windows:

- A C function that calls any libc function `XSUB.h` redefines (`malloc`,
  `free`, `realloc`, `calloc`, `fopen`, `getenv`, `exit`, sockets, …) takes
  `pTHX_`, even a static helper that touches no SV: every Windows perl has
  `PERL_IMPLICIT_SYS`, which routes them through `my_perl`. 0.3414's
  `alloc_probe_ok()` lacked it and failed on Strawberry. `./test.all.perls.pl`
  emulates this with a `win32-sys` row (a syntax check of `LikeR.c` with
  `-DPERL_IMPLICIT_SYS -DSOCKET=int` against the newest threaded perl); a FAIL
  there is a real Windows build failure. Do not pass `--no-win-sys` for a
  change to `LikeR.xs`.
- No `<unistd.h>`, `<sys/time.h>`, `<strings.h>`, `strcasecmp`/`strncasecmp`
  (use `str_ieq_ascii()`), `fork`, or hardcoded `/` in paths.
- Printf lengths are perl's: `%" UVuf "`, `%" IVdf "`, `%" NVgf "`,
  `%" UVxf "`; never `%zu`, `%zd`, `%lld`, `%llu` — cast to `UV`/`IV`.
- A test never builds a shell command line: backticks, `qx//` and one-argument
  `system()` go through `cmd.exe`, which does not treat `'` as a quote. Use
  list-form `system()`, put the program in a file, and pass paths in `@ARGV`.
  `t/arg.crash.regressions.t` is the worked example (learned in 0.317).
- A test never compares a byte count on disk (`-s`, `stat`) with the length of
  text that went through a text-mode handle: `:crlf` writes each `\n` as
  `\r\n` and reads it back as one byte, so the two differ by one per line.
  Read the file back through the same layers and compare the text, or open
  both ends `:raw`. Chem::Structure::Parser 0.040 failed on every Strawberry
  smoker this way. `-s` as a non-empty check, or against the same file's own
  earlier size, is fine.
- Never depend on renaming over, or deleting, a file that someone may hold
  open: Windows refuses while any handle has it open (the C runtime never
  asks for `FILE_SHARE_DELETE`), and a `File::Temp` object, or a
  `tempfile()` made with `UNLINK`, holds one until `END`. Opening the file to
  write in place is allowed, so keep that as the fallback, as
  `wt_target_copy()` does for `write_table()`. 0.3428's atomic replace had
  none, and failed three test files on every Strawberry smoker. On Linux, an
  `LD_PRELOAD` `rename()` that returns EACCES while the process has the
  target open (in `/proc/self/fd`) reproduces it.

Solaris, illumos and the BSDs:

- Add to the `_GNU_SOURCE` / `__EXTENSIONS__` block at the top of
  `LikeR.xs`, never around it; `__EXTENSIONS__` exposes POSIX and XPG
  symbols on Solaris under strict C99.
- No glibc- or BSD-only libc: `qsort_r` (argument order differs; use
  `sortsv`/`sortsv_flags` or `qsort` with a file-static comparator),
  `memmem`, `strchrnul`, `strndup`, `asprintf`, `getline`, `fmemopen`,
  `reallocarray`, `explicit_bzero`, `arc4random`, `exp10`, `sincos`.
- Stay within C99 `<math.h>` (`lgamma`, `log2`, `round`, `trunc` and `fmax`
  are all there); anything beyond (`lgamma_r`, `j0`, `tgammal`,
  `lgammaq`) needs a guard and a fallback. OpenBSD is thin on long-double
  variants.
- A randomness path keeps its documented fallback for when `/dev/urandom`
  cannot be opened.
- Assume nothing about struct field order, `char` signedness or alignment;
  never cast a `char *` to a wider pointer type and dereference it.

### Back-compatible to perl 5.10

`dist.ini` declares `perl = 5.010` and `lib/Stats/LikeR.pm` says
`require 5.010`.

- Perl code (module and tests) uses 5.10 syntax only: no signatures, postfix
  dereference, unimported `say`, `isa` or chained comparisons. Tests use the
  existing header (`require 5.010; use strict; use warnings FATAL => 'all';`),
  not `use 5.044`. The unshipped root helpers (`test.all.perls.pl`,
  `xs.check.pl`) are exempt.
- Perl C API newer than 5.10 comes through `ppport.h`, regenerated when a new
  one is used; if `ppport.h` lacks it, write the fallback.
- Do not use API that was merely renamed later, or rely on behavior changed
  after 5.10 (hash order, `SvPV` return guarantees, COW semantics) without
  checking on `perl-5.10.1`.
- Raising `MIN_PERL_VERSION` is the maintainer's decision, never a way to fix
  a build error. If a change truly needs a newer perl, say so and stop.

#### Never evaluate a `qr//` inside a leak-test block

On 5.10.0, which smokers still run and nothing here reproduces, `pp_qr()`
leaks one SV (`PV "Regexp"`) per `qr//` evaluated; 5.10.1 fixed it. 0.315
failed `t/cfilter.t` and `t/filter_match.t` there. Compile the pattern into a
lexical *before* the `no_leaks_ok` block, pass that in, and say why in a
comment; `cfilter` and `col()->match` take a precompiled `qr//` unchanged.
