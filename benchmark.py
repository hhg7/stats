#!/usr/bin/env python3
import time
import tracemalloc
import csv
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn import metrics

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
    'cor': ('df.corr', lambda: df[['x', 'y']].corr()),
    'cov': ('df.cov', lambda: df[['x', 'y']].cov()),
    
    # Distributions & Tests
    'rnorm': ('np.random.randn', lambda: np.random.randn(1000)),
    'runif': ('np.random.rand', lambda: np.random.rand(1000)),
    't_test': ('scipy.stats.ttest_ind', lambda: stats.ttest_ind(df[df['cat2']=='X']['x'], df[df['cat2']=='Y']['x'])),
    'wilcox_test': ('scipy.stats.mannwhitneyu', lambda: stats.mannwhitneyu(df[df['cat2']=='X']['x'], df[df['cat2']=='Y']['x'])),
    'chisq_test': ('scipy.stats.chi2_contingency', lambda: stats.chi2_contingency(pd.crosstab(df['cat1'], df['cat2']))),
    'shapiro_test': ('scipy.stats.shapiro', lambda: stats.shapiro(df['x'][:5000])), # shapiro has size limits
    'binom_test': ('scipy.stats.binomtest', lambda: stats.binomtest(500, 1000, 0.5)),
    
    # Data Manipulation (Pandas)
    'filter': ('df.query', lambda: df.query("x > 0")),
    'select_cols': ('df.__getitem__', lambda: df[['x', 'cat1']]),
    'drop_cols': ('df.drop', lambda: df.drop(columns=['y'])),
    'rename_cols': ('df.rename', lambda: df.rename(columns={'cat1': 'Category_1'})),
    'dropna': ('df.dropna', lambda: df_missing.dropna()),
    'fillna': ('df.fillna', lambda: df_missing.fillna(0)),
    'drop_duplicates': ('df.drop_duplicates', lambda: df.drop_duplicates()),
    'group_by': ('df.groupby', lambda: df.groupby('cat1').mean(numeric_only=True)),
    'concat': ('pd.concat', lambda: pd.concat([df, df])),
    'merge': ('pd.merge', lambda: pd.merge(df, df, left_index=True, right_index=True, suffixes=('_1', '_2'))),
    'value_counts': ('df.value_counts', lambda: df['cat1'].value_counts()),
    'pivot_table': ('pd.pivot_table', lambda: pd.pivot_table(df, values='x', index='cat1', columns='cat2', aggfunc='mean')),
    
    # Modeling & Metrics
    'lm': ('statsmodels.OLS', lambda: smf.ols('y ~ x + C(cat1)', data=df).fit()),
    'glm': ('statsmodels.GLM', lambda: smf.logit('binary ~ x + y', data=df).fit()),
    'auc': ('sklearn.metrics.auc', lambda: metrics.auc(*metrics.roc_curve(y_true, y_scores)[:2])),
    'prcomp': ('sklearn.decomposition.PCA', lambda: __import__('sklearn.decomposition').decomposition.PCA(n_components=2).fit(df[['x', 'y']]))
}

# 3. Execution Engine
results = []
runs = 7

print("Running Python benchmarks...")
for liker_func, (py_sub, func) in benchmarks.items():
    for run in range(runs):
        # Force garbage collection before tracking
        import gc
        gc.collect()
        
        tracemalloc.start()
        start_time = time.perf_counter()
        
        try:
            _ = func()
        except Exception as e:
            print(f"Error in {liker_func}: {e}")
            
        end_time = time.perf_counter()
        current, peak_mem = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        elapsed_time = end_time - start_time
        
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