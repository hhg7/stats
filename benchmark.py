#!/usr/bin/env python3
# The Python side of benchmark.R and benchmark.pl: the same 32 operations, on
# the same 10,000-row mock frame, seven runs each, written out as
# python_benchmarks.tsv with the same four columns the other two scripts write.
#
#     python3 benchmark.py
#
# Where the mapping is not one-to-one, a comment on the entry says why.  The
# rule everywhere is that the three languages must be asked for the *same
# result* by each one's idiomatic route -- not for whichever call happens to
# look most like the others' spelling.  Two ways that gets broken, both of
# which this file has had to correct for:
#
#   * asking pandas for more than the others.  df[['x','y']].corr() returns a
#     2x2 matrix where Stats::LikeR's cor() returns one number, and
#     groupby(...).mean() averages every numeric column where R's summarise()
#     averages one.  Both were doing 1.6-2.8x the work of their counterparts.
#   * asking pandas for less.  Merging on the index is not the join the other
#     two perform; they join on a key column, which is about twice the work.
import gc
import time
import tracemalloc
import csv
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn import metrics
from sklearn.decomposition import PCA

# 1. Setup Mock Data
np.random.seed(42)
n = 10000
df = pd.DataFrame({
    'x': np.random.randn(n),
    'y': np.random.randn(n) * 2 + 5,
    'cat1': np.random.choice(['A', 'B', 'C'], n),
    'cat2': np.random.choice(['X', 'Y'], n),
    'binary': np.random.binomial(1, 0.5, n)
})
df_missing = df.copy()
df_missing.loc[10:20, 'x'] = np.nan

# merge needs a key column: Stats::LikeR has no index to join on, so all three
# scripts join two copies of the frame on an explicit id.
df_id = df.copy()
df_id['id'] = np.arange(1, n + 1)

y_true = np.random.binomial(1, 0.5, n)
y_scores = np.random.rand(n)

# 2. Define Function Mappings (Stats::LikeR -> Python equivalent)
# We use lambdas to pass the pre-generated data into the functions
benchmarks = {
    # Basic Stats
    'mean': ('np.mean', lambda: np.mean(df['x'])),
    'median': ('np.median', lambda: np.median(df['x'])),
    'sd': ('np.std', lambda: np.std(df['x'], ddof=1)),
    'var': ('np.var', lambda: np.var(df['x'], ddof=1)),
    'min': ('np.min', lambda: np.min(df['x'])),
    'max': ('np.max', lambda: np.max(df['x'])),
    'quantile': ('np.quantile', lambda: np.quantile(df['x'], [0.25, 0.5, 0.75])),
    # one coefficient, not a 2x2 matrix: Stats::LikeR's cor/cov take two vectors
    'cor': ('pd.Series.corr', lambda: df['x'].corr(df['y'])),
    'cov': ('np.cov', lambda: np.cov(df['x'], df['y'])[0, 1]),

    # Distributions & Tests
    'rnorm': ('np.random.randn', lambda: np.random.randn(1000)),
    'runif': ('np.random.rand', lambda: np.random.rand(1000)),
    # Welch, which is what R's t.test and Stats::LikeR's t_test both default to
    't_test': ('scipy.stats.ttest_ind', lambda: stats.ttest_ind(df[df['cat2']=='X']['x'], df[df['cat2']=='Y']['x'], equal_var=False)),
    'wilcox_test': ('scipy.stats.mannwhitneyu', lambda: stats.mannwhitneyu(df[df['cat2']=='X']['x'], df[df['cat2']=='Y']['x'])),
    # building the table is inside the timed call, as it is for the other two
    'chisq_test': ('scipy.stats.chi2_contingency', lambda: stats.chi2_contingency(pd.crosstab(df['cat1'], df['cat2']))),
    'shapiro_test': ('scipy.stats.shapiro', lambda: stats.shapiro(df['x'][:5000])), # shapiro has size limits
    'binom_test': ('scipy.stats.binomtest', lambda: stats.binomtest(500, 1000, 0.5)),

    # Data Manipulation (Pandas)
    # boolean mask, the pandas idiom; .query() sends the string through the
    # expression parser first, which is ~2x the cost and has no counterpart in
    # filter($df, col('x') > 0) -- that predicate is compiled once, not parsed
    # per call
    'filter': ('df.__getitem__ (mask)', lambda: df[df['x'] > 0]),
    'select_cols': ('df.__getitem__', lambda: df[['x', 'cat1']]),
    'drop_cols': ('df.drop', lambda: df.drop(columns=['y'])),
    'rename_cols': ('df.rename', lambda: df.rename(columns={'cat1': 'Category_1'})),
    'dropna': ('df.dropna', lambda: df_missing.dropna()),
    'fillna': ('df.fillna', lambda: df_missing.fillna(0)),
    'drop_duplicates': ('df.drop_duplicates', lambda: df.drop_duplicates()),
    # one column's mean per group, as in R's summarise(mean_x = mean(x)) and in
    # the Perl group_by + mean; .mean(numeric_only=True) would average three
    'group_by': ('df.groupby', lambda: df.groupby('cat1')['x'].mean()),
    'concat': ('pd.concat', lambda: pd.concat([df, df])),
    'merge': ('pd.merge', lambda: pd.merge(df_id, df_id, on='id', suffixes=('_1', '_2'))),
    'value_counts': ('df.value_counts', lambda: df['cat1'].value_counts()),
    'pivot_table': ('pd.pivot_table', lambda: pd.pivot_table(df, values='x', index='cat1', columns='cat2', aggfunc='mean')),

    # Modeling & Metrics
    'lm': ('statsmodels.OLS', lambda: smf.ols('y ~ x + C(cat1)', data=df).fit()),
    # IRLS, the algorithm behind R's glm(family=binomial) and Stats::LikeR's
    # glm(family => 'binomial'); smf.logit would solve it by Newton-Raphson
    'glm': ('statsmodels.GLM', lambda: smf.glm('binary ~ x + y', data=df, family=sm.families.Binomial()).fit()),
    'auc': ('sklearn.metrics.auc', lambda: metrics.auc(*metrics.roc_curve(y_true, y_scores)[:2])),
    # fit_transform, because R's prcomp and Stats::LikeR's both return the
    # scores; a bare .fit() stops one matrix multiply short of them
    'prcomp': ('sklearn.decomposition.PCA', lambda: PCA(n_components=2).fit_transform(df[['x', 'y']]))
}

# 3. Execution Engine
#
# Time and memory are measured in SEPARATE runs.  tracemalloc hooks every
# allocation, so a call made while it is running takes 1.4x to 3x as long as
# the same call made without it -- for the statsmodels fits, three times.
# benchmark.R takes its gc() readings either side of Sys.time(), and
# benchmark.pl reads /proc either side of the timer in a forked child; neither
# has a profiler attached while the clock is running, and neither can this one,
# or the column measures tracemalloc rather than pandas.  The cost is that each
# operation is executed 2 * runs times instead of runs.
#
# What the RAM number contains: every allocation Python's allocator sees during
# the call, at its high-water mark, including NumPy and pandas buffers (NumPy
# routes them through PyDataMem, which tracemalloc traces).  It does not see
# plain malloc from BLAS or SciPy's Fortran kernels.  Like the gc() figure in
# benchmark.R and the RSS delta in benchmark.pl, it is the right order of
# magnitude rather than an exact ledger, and the three are not the same
# quantity -- compare them across operations, not to each other.
results = []
runs = 7

print("Running Python benchmarks...")
for liker_func, (py_sub, func) in benchmarks.items():
    for run in range(runs):
        # --- timed run: nothing else attached ---
        gc.collect()
        start_time = time.perf_counter()
        try:
            _ = func()
            elapsed_time = time.perf_counter() - start_time
        except Exception as e:
            elapsed_time = float('nan')
            print(f"Error in {liker_func}: {e}")

        # --- traced run: same call, no clock ---
        gc.collect()
        tracemalloc.start()
        try:
            _ = func()
            _, peak_mem = tracemalloc.get_traced_memory()
        except Exception:
            peak_mem = float('nan')
        finally:
            tracemalloc.stop()

        results.append({
            'Stats::LikeR function': liker_func,
            'Python sub': py_sub,
            'time': elapsed_time,
            'RAM': peak_mem
        })

# 4. Output
output_file = 'python_benchmarks.tsv'
with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['Stats::LikeR function', 'Python sub', 'time', 'RAM'], delimiter='\t')
    writer.writeheader()
    writer.writerows(results)

print(f"Done. Python results written to {output_file}")
