#!/usr/bin/env python3
"""Time read/write of t/HepatitisCdata.csv over 10 runs each.

read_table  -> pandas.read_csv   (csv-parsing form; his read_table -> _parse_csv_file)
write_table -> DataFrame.to_csv

Reads t/HepatitisCdata.csv; writes the same data to a temp file in /tmp.
Emits python.read_table.tsv and python.write_table.tsv, one row per run.
Optional argv[1] overrides the input path.
"""
import os
import sys
import tempfile
import time

import pandas as pd


def save_times(fname, times):
    with open(fname, 'w') as fh:
        fh.write('run\tseconds\n')
        for i, t in enumerate(times):
            fh.write(f'{i}\t{t:.9f}\n')


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else 't/HepatitisCdata.csv'
    d = pd.read_csv(src)

    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.csv', dir='/tmp', delete=False)
    tmp.close()
    path = tmp.name

    write = []
    for _ in range(10):
        t0 = time.perf_counter()
        d.to_csv(path, index=False)
        t1 = time.perf_counter()
        write.append(t1 - t0)
    save_times('python.write_table.tsv', write)

    read = []
    for _ in range(10):
        t0 = time.perf_counter()
        _ = pd.read_csv(src)
        t1 = time.perf_counter()
        read.append(t1 - t0)
    save_times('python.read_table.tsv', read)

    os.unlink(path)


if __name__ == '__main__':
    main()
