#!/usr/bin/env python3
"""The Python side of plot.scaling.pl: the same functions, the same size ladder.

Where benchmark.py asks "how long does each function take on one 10,000-row
frame", this asks "what shape is that curve".  Each call is timed at a ladder
of sizes half a decade apart, seven times at each size, and written out for
plot.scaling.pl to draw against the same measurements from Perl and R.  On a
log-log axis the slope of the line is the exponent.

    perl plot.scaling.pl --data                 # once, writes the fixtures
    perl -Iblib/arch -Iblib/lib \
         plot.scaling.pl --measure              # -> perl_scaling.tsv
    python3 scale.py                            # -> python_scaling.tsv
    Rscript scale.R                             # -> r_scaling.tsv
    perl plot.scaling.pl --plot                 # -> scaling.*.svg

Environment:

    SCALE_DIR    where --data put the fixtures (/tmp/likeR.scaling)
    SCALE_RUNS   runs per (function, size); default 7
    SCALE_CAP    seconds; once one run of a function takes longer than this,
                 that function is not tried at any larger size.  Default 4.
    SCALE_MAX_N  hard ceiling on the row count, for a quick partial run
    SCALE_TARGET seconds a single measurement should span; a call faster than
                 this is repeated until it does.  Default 0.002.
    SCALE_CPU    the CPU this run pins itself to; default 0, "none" to let the
                 scheduler place it.

On a hybrid CPU -- Intel's P-core/E-core parts, and the big.LITTLE ARM designs
-- a process landing on an E-core reads 1.5 to 2 times slower than the same
loop on a P-core, which is a wider band than most of the differences these
plots are drawn to show.  So this pins itself to one CPU rather than asking to
be started under "taskset"; see pin_to_one_cpu() below.  plot.scaling.pl
--measure and scale.R do the same, and all three default to CPU 0, so the three
languages land on the same core without anyone having to remember.  SCALE_CPU
moves it, and has to be given the same value in all three.

The rule from benchmark.py carries over unchanged: each language is asked for
the *same result* by its own idiomatic route, not for whichever call happens to
be spelled most like the others.  Two ways that gets broken here, both guarded
against below:

  * asking pandas for less.  numpy's .T is a view, not a transposition, and
    df[['x','y']] is a view onto the parent's blocks; Stats::LikeR's transpose
    and select_cols both build a new structure, so the .copy() is what makes
    the two calls the same job.
  * asking pandas for more.  DataFrame.corr() would return a 2x2 matrix where
    cor() returns one number, and groupby().mean() would average every numeric
    column where group_by + mean averages one.

The 'function' strings must match plot.scaling.pl's exactly -- that is the key the
plot joins the three files on.  A name here with no counterpart there simply
draws one fewer line in that panel.
"""
import csv
import os
import sys
import time

import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn import metrics

DIR = os.environ.get('SCALE_DIR', '/tmp/likeR.scaling')
RUNS = int(os.environ.get('SCALE_RUNS', 7))
CAP = float(os.environ.get('SCALE_CAP', 4))
MAX_N = int(os.environ.get('SCALE_MAX_N', 0))
TARGET = float(os.environ.get('SCALE_TARGET', 0.002))
MAX_REPS = 10_000

# Size ladders.  Half-decade steps; plot.scaling.pl and scale.R carry the same
# three lists, and the fixture row counts are IO_N.  Change one, change all three.
VEC_N = [1_000, 3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000]
IO_N = [1_000, 3_000, 10_000, 30_000, 100_000, 300_000]
FRAME_N = [1_000, 3_000, 10_000, 30_000, 100_000, 300_000]

if MAX_N:
    VEC_N = [n for n in VEC_N if n <= MAX_N]
    IO_N = [n for n in IO_N if n <= MAX_N]
    FRAME_N = [n for n in FRAME_N if n <= MAX_N]


# ---------------------------------------------------------------------------
# Input builders -- one per figure, called once per size
# ---------------------------------------------------------------------------
def build_vector(n):
    rng = np.random.default_rng(42)
    return {
        'x': rng.standard_normal(n),
        'y': rng.standard_normal(n) * 2 + 5,
        'label': rng.binomial(1, 0.5, n),
        'n': n,
    }


def build_io(n):
    f = {
        'num_csv': os.path.join(DIR, 'num.%d.csv' % n),
        'mix_csv': os.path.join(DIR, 'mix.%d.csv' % n),
        'mix_tsv': os.path.join(DIR, 'mix.%d.tsv' % n),
        'mix_xlsx': os.path.join(DIR, 'mix.%d.xlsx' % n),
        'out': os.path.join(DIR, 'out.python.%d.tmp' % os.getpid()),
        # to_excel picks its engine off the extension, as write_table does
        'out_xlsx': os.path.join(DIR, 'out.python.%d.xlsx' % os.getpid()),
    }
    for k in ('num_csv', 'mix_csv', 'mix_tsv', 'mix_xlsx'):
        if not os.path.isfile(f[k]):
            sys.exit('missing fixture "%s"; run "perl plot.scaling.pl --data" first'
                     % f[k])
    # The frame handed to to_csv is read in, not synthesized, so all three
    # languages write out the same table.
    f['df'] = pd.read_csv(f['mix_csv'])
    return f


def build_frame(n):
    rng = np.random.default_rng(42)
    df = pd.DataFrame({
        'x': rng.standard_normal(n),
        'y': rng.standard_normal(n) * 2 + 5,
        'cat1': rng.choice(['A', 'B', 'C'], n),
        'cat2': rng.choice(['X', 'Y'], n),
        'binary': rng.binomial(1, 0.5, n),
    })
    # merge needs a key column: a Stats::LikeR frame has no index to join on,
    # so all three scripts join on an explicit id, as benchmark.py does.
    df_id = df.copy()
    df_id['id'] = np.arange(1, n + 1)
    return {
        'df': df,
        'df_id': df_id,
        # aoh2hoa's input: one record per row, which is what
        # pd.DataFrame(records) is handed.
        'records': df.to_dict('records'),
        'arr': df[['x', 'y', 'binary']].to_numpy(),
        'n': n,
    }


BUILD = {'vector': build_vector, 'io': build_io, 'frame': build_frame}

# ---------------------------------------------------------------------------
# The benchmarks: (figure, function, call, body)
# ---------------------------------------------------------------------------
# 'function' is the join key with plot.scaling.pl and scale.R; 'call' is only recorded
# for the reader.
BENCHMARKS = [
    # --- reductions over one numeric vector --------------------------------
    ('vector', 'sum', 'np.sum(x)', lambda d: np.sum(d['x'])),
    ('vector', 'min', 'np.min(x)', lambda d: np.min(d['x'])),
    ('vector', 'max', 'np.max(x)', lambda d: np.max(d['x'])),
    ('vector', 'mean', 'np.mean(x)', lambda d: np.mean(d['x'])),
    ('vector', 'median', 'np.median(x)', lambda d: np.median(d['x'])),
    ('vector', 'sd', 'np.std(x, ddof=1)', lambda d: np.std(d['x'], ddof=1)),
    ('vector', 'var', 'np.var(x, ddof=1)', lambda d: np.var(d['x'], ddof=1)),
    ('vector', 'quantile', 'np.quantile(x, [.25,.5,.75])',
     lambda d: np.quantile(d['x'], [0.25, 0.5, 0.75])),
    # one number, not the 2x2 matrix np.corrcoef/df.corr would build
    ('vector', 'cor', 'np.corrcoef(x, y)[0, 1]',
     lambda d: np.corrcoef(d['x'], d['y'])[0, 1]),
    ('vector', 'cov', 'np.cov(x, y)[0, 1]',
     lambda d: np.cov(d['x'], d['y'])[0, 1]),
    ('vector', 'skew', 'scipy.stats.skew(x)', lambda d: stats.skew(d['x'])),
    ('vector', 'kurtosis', 'scipy.stats.kurtosis(x)',
     lambda d: stats.kurtosis(d['x'])),

    # --- transforms that return something the size of their input ----------
    ('transform', 'rank', 'scipy.stats.rankdata(x)',
     lambda d: stats.rankdata(d['x'])),
    ('transform', 'uniq', 'pd.unique(x)', lambda d: pd.unique(d['x'])),
    ('transform', 'scale', 'scipy.stats.zscore(x, ddof=1)',
     lambda d: stats.zscore(d['x'], ddof=1)),
    ('transform', 'sample', 'rng.choice(x, n // 10, replace=False)',
     lambda d: np.random.default_rng(1).choice(d['x'], d['n'] // 10 + 1,
                                               replace=False)),
    ('transform', 'seq', 'np.arange(1, n + 1)',
     lambda d: np.arange(1, d['n'] + 1)),
    ('transform', 'auc', 'sklearn.metrics.roc_auc_score(label, y)',
     lambda d: metrics.roc_auc_score(d['label'], d['y'])),

    # --- read_table and write_table, over four inputs each -----------------
    # pandas is columnar, so read_csv is the honest counterpart of
    # read_table(output.type => 'hoa') and appears in that panel too; the
    # row-record panels get csv.DictReader, which is what actually produces the
    # shape read_table returns by default.
    ('io', 'read_table (csv, numeric)', 'pd.read_csv(num.csv)',
     lambda d: pd.read_csv(d['num_csv'])),
    ('io', 'read_table (csv, mixed)', 'list(csv.DictReader(mix.csv))',
     lambda d: _dictread(d['mix_csv'], ',')),
    ('io', 'read_table (tsv, mixed)', 'list(csv.DictReader(mix.tsv, "\\t"))',
     lambda d: _dictread(d['mix_tsv'], '\t')),
    ('io', 'read_table (csv, hoa)', 'pd.read_csv(mix.csv)',
     lambda d: pd.read_csv(d['mix_csv'])),
    # The same table as the three panels above, out of a spreadsheet.
    # read_excel defaults to openpyxl for .xlsx, which is what a pandas user
    # gets without asking for anything; it is a pure-python parser and reads an
    # order of magnitude slower than the csv path, so SCALE_CAP will usually
    # stop this panel a size or two before the others.
    ('io', 'read_table (xlsx)', 'pd.read_excel(mix.xlsx)',
     lambda d: pd.read_excel(d['mix_xlsx'])),
    ('io', 'write_table (csv, hoa)', 'df.to_csv(f, index=False)',
     lambda d: d['df'].to_csv(d['out'], index=False)),
    ('io', 'write_table (tsv, hoa)', 'df.to_csv(f, sep="\\t", index=False)',
     lambda d: d['df'].to_csv(d['out'], sep='\t', index=False)),
    ('io', 'write_table (csv, aoa)', 'df.to_csv(f, index=False)',
     lambda d: d['df'].to_csv(d['out'], index=False)),
    ('io', 'write_table (csv, row.names)', 'df.to_csv(f, index=True)',
     lambda d: d['df'].to_csv(d['out'], index=True)),
    # openpyxl again, on the way out.  R is absent from this panel: it has no
    # xlsx writer installed here and base R has none at all.
    ('io', 'write_table (xlsx, hoa)', 'df.to_excel(f, index=False)',
     lambda d: d['df'].to_excel(d['out_xlsx'], index=False)),

    # --- whole-frame operations --------------------------------------------
    ('frame', 'filter', 'df[df.x > 0]', lambda d: d['df'][d['df']['x'] > 0]),
    # No .copy() here, unlike transpose below.  select_cols on a HoA returns a
    # new hash holding the *same* array references -- the columns are aliased,
    # not copied, which is why its panel is a flat line -- so the pandas view
    # is the matching job and .copy() would be asking pandas for more.  R's
    # subset does materialize, and that difference is what the panel shows.
    ('frame', 'select_cols', "df[['x','cat1']]",
     lambda d: d['df'][['x', 'cat1']]),
    # one column averaged per group, not every numeric one
    ('frame', 'group_by + mean', "df.groupby('cat1')['x'].mean()",
     lambda d: d['df'].groupby('cat1')['x'].mean()),
    ('frame', 'merge', "df_id.merge(df_id, on='id', how='inner')",
     lambda d: d['df_id'].merge(d['df_id'], on='id', how='inner')),
    ('frame', 'value_counts', "df['cat1'].value_counts()",
     lambda d: d['df']['cat1'].value_counts()),
    ('frame', 'drop_duplicates', 'df.drop_duplicates()',
     lambda d: d['df'].drop_duplicates()),
    # .copy(): numpy's .T is a view, transpose() builds a new AoA
    ('frame', 'transpose', 'arr.T.copy()', lambda d: d['arr'].T.copy()),
    ('frame', 'aoh2hoa', 'pd.DataFrame(records)',
     lambda d: pd.DataFrame(d['records'])),
]


def _dictread(path, sep):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter=sep))


LADDER = {'vector': VEC_N, 'transform': VEC_N, 'io': IO_N, 'frame': FRAME_N}
BUILDER_FOR = {'vector': 'vector', 'transform': 'vector', 'io': 'io',
               'frame': 'frame'}


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------
# Timed in-process, unlike plot.scaling.pl, which forks: CPython's allocator does not
# hand a later run a materially different starting point, and the numpy and
# pandas import cost that a fork per run would repeat dwarfs everything being
# measured.  What does need controlling is the first call -- pandas resolves
# dtypes and compiles paths lazily -- so every function gets one untimed warm-up
# at every size before its timed runs, which also warms the page cache for the
# file panels.
#
# A call faster than TARGET is repeated until the pair of clock readings spans
# TARGET and the total divided by the count.  perf_counter resolves
# nanoseconds, so this is not about the clock the way it is in scale.R; it is
# about the interpreter, whose dispatch overhead around a two-microsecond
# np.sum is a large fraction of what would otherwise be recorded.  The repeat
# count is chosen from a second untimed call, so that all three scripts average
# over the same span of work.
# The CPU this run is pinned to, fixed before the first measurement.
#
# Saying "run this under taskset" in a comment is not the same as running it
# under taskset.  An unpinned rerun of an identical build read Stats::LikeR's
# min() at 2.75 ms against 1.43 pinned, and mean() at 2.61 against 1.35 -- a
# factor of 1.9 on a 13th-generation Core i7, landing on whichever functions
# happened to draw an E-core rather than on all of them evenly.  A reading that
# moves by 1.9x depending on where the scheduler put the process is not a
# measurement, and a cross-language plot drawn from three of them unpinned is
# comparing schedulers.
#
# os.sched_setaffinity is the whole mechanism: it is Linux-only, but it is in
# the standard library, so unlike plot.scaling.pl and scale.R this does not have
# to shell out to taskset(1).
#
# Every failure is survivable and none of them stop the run:
#
#   * No sched_getaffinity -- macOS, Windows -- means the affinity cannot be
#     read, let alone set, so this says so once and carries on.
#   * A setaffinity that is refused (a cgroup or a taskset that already excluded
#     the CPU asked for) leaves the run on the CPUs it had.
#   * A process already pinned to one CPU is left on it, so
#     "taskset -c 5 python3 scale.py" keeps CPU 5 and this does not overrule
#     the operator.
#
# CPU 0 is the default because the P-cores are enumerated first on the hybrid
# parts, so it is a fast core where the choice matters at all, and because it is
# what the other two scripts default to.
def pin_to_one_cpu():
    want = os.environ.get('SCALE_CPU', '0')
    if want in ('none', ''):
        print('SCALE_CPU=none: leaving the scheduler to place the measurements')
        return
    try:
        cpu = int(want)
    except ValueError:
        cpu = -1
    if cpu < 0:
        sys.exit("SCALE_CPU must be a CPU number or 'none', not %r" % want)

    try:
        before = os.sched_getaffinity(0)
    except AttributeError:
        print('this platform cannot report a CPU affinity, so the measurements '
              'are not pinned; on a hybrid CPU, start this under the '
              "platform's affinity tool")
        return
    if len(before) == 1:
        print('already pinned to CPU %d; keeping it' % next(iter(before)))
        return
    try:
        os.sched_setaffinity(0, {cpu})
    except OSError as exc:
        print('could not pin to CPU %d, staying on the %d it had: %s'
              % (cpu, len(before), exc), file=sys.stderr)
        return
    print('pinned to CPU %d (was %d CPUs); plot.scaling.pl --measure and '
          'scale.R must be on the same one' % (cpu, len(before)))


def measure(body, data):
    """One reading: seconds per call, averaged over however many calls fit."""
    t0 = time.perf_counter()
    body(data)
    one = time.perf_counter() - t0
    reps = min(MAX_REPS, int(TARGET / one) + 1) if one > 0 else MAX_REPS
    t0 = time.perf_counter()
    for _ in range(reps):
        body(data)
    return (time.perf_counter() - t0) / reps, reps


def main():
    pin_to_one_cpu()
    results = []
    too_slow = set()

    by_figure = {}
    for figure, name, call, body in BENCHMARKS:
        by_figure.setdefault(figure, []).append((name, call, body))

    for figure in ('vector', 'transform', 'io', 'frame'):
        if figure not in by_figure:
            continue
        for n in LADDER[figure]:
            todo = [b for b in by_figure[figure] if b[0] not in too_slow]
            if not todo:
                continue
            data = BUILD[BUILDER_FOR[figure]](n)
            for name, call, body in todo:
                try:
                    body(data)                    # warm-up, not recorded
                except Exception as exc:          # noqa: BLE001 - report and move on
                    print('%s at n=%d: %s' % (name, n, exc), file=sys.stderr)
                    too_slow.add(name)
                    continue

                slowest = 0.0
                reps = 0
                for run in range(RUNS):
                    elapsed, reps = measure(body, data)
                    slowest = max(slowest, elapsed)
                    results.append((figure, name, call, n, run, elapsed))
                print('%-9s %-30s n=%-8d %.6f s%s'
                      % (figure, name, n, slowest,
                         ' (x%d)' % reps if reps > 1 else ''))
                if slowest > CAP:
                    too_slow.add(name)
            del data

    # A curve that ends early is not a curve that was never measured, and the
    # plot cannot tell you which one you are looking at.
    if too_slow:
        print('Stopped early (a run exceeded %g s, or it failed): %s'
              % (CAP, ', '.join(sorted(too_slow))))

    for out in (os.path.join(DIR, 'out.python.%d.tmp' % os.getpid()),
                os.path.join(DIR, 'out.python.%d.xlsx' % os.getpid())):
        if os.path.isfile(out):
            os.unlink(out)

    with open('python_scaling.tsv', 'w') as fh:
        fh.write('figure\tfunction\tcall\tn\trun\tseconds\n')
        for row in results:
            fh.write('%s\t%s\t%s\t%d\t%d\t%.9f\n' % row)
    print('Done. %d measurements written to python_scaling.tsv' % len(results))


if __name__ == '__main__':
    main()
