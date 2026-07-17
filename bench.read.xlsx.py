#!/usr/bin/env python3
"""Time reading an .xlsx file over 10 runs and report memory usage.

Equivalent of benchmark.pl's timing loop, but for read_table on .xlsx:
  read xlsx -> pandas.read_excel (pandas' default engine)

Reads titanic.more.complete.xlsx by default (argv[1] overrides).
Emits xlsx.read.Python.tsv, one row per run (run\tseconds).
Prints memory usage to stdout:
  - peak Python allocations during a read (tracemalloc)
  - in-memory DataFrame size (deep)
  - process peak resident set size (ru_maxrss)
"""
import resource
import sys
import time
import tracemalloc

import pandas as pd


def save_times(fname, times):
    with open(fname, 'w') as fh:
        fh.write('run\tseconds\n')
        for i, t in enumerate(times):
            fh.write(f'{i}\t{t:.9f}\n')


def human(nbytes):
    for unit in ('B', 'KiB', 'MiB', 'GiB'):
        if nbytes < 1024 or unit == 'GiB':
            return f'{nbytes:.2f} {unit}'
        nbytes /= 1024


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else 'titanic.more.complete.xlsx'

    read = []
    for _ in range(10):
        t0 = time.perf_counter()
        d = pd.read_excel(src)
        t1 = time.perf_counter()
        read.append(t1 - t0)
    save_times('xlsx.read.Python.tsv', read)

    # peak Python allocations for a single read
    tracemalloc.start()
    d = pd.read_excel(src)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    df_bytes = int(d.memory_usage(deep=True).sum())
    max_rss_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss  # KiB on Linux

    print(f'file                : {src}')
    print(f'shape               : {d.shape[0]} rows x {d.shape[1]} cols')
    print(f'median read time    : {sorted(read)[len(read) // 2]:.6f} s')
    print(f'DataFrame in-memory : {human(df_bytes)}')
    print(f'peak Python alloc   : {human(peak)}')
    print(f'process peak RSS    : {human(max_rss_kb * 1024)}')


if __name__ == '__main__':
    main()
