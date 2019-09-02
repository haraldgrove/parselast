#! /usr/bin/env python

"""
Calculate the coverage at each position for the target sequence
Only recommended for targets less than 10kb.
"""

import sys
import pandas as pd


def roundup(x, n):
    N = pow(10, n)
    return x if x % N == 0 else x + N - x % N


def calc_coverage(filename):
    """
    """
    print("Reading {}".format(filename))
    df = pd.read_table(
        filename, comment="#", header=0, delim_whitespace=True, error_bad_lines=False
    )
    cov = {}
    for (index, row) in df.iterrows():
        if row['idpct'] < 100 or row['covpct'] < 50:
            continue
        if row["name1"] not in cov:
            cov[row["name1"]] = [0] * row['seqSize1']
        for pos in range(row['start1'], row['end1']):
            cov[row['name1']][pos] += 1
    outputname = '{}.coverage.txt'.format(filename.rsplit('.', 1)[0])
    with open(outputname, 'w') as fout:
        for target, coverage in cov.items():
            for pos, c in enumerate(coverage):
                fout.write('{}\t{}\t{}\n'.format(target, pos, c))


def main():
    for filename in sys.argv[1:]:
        try:
            calc_coverage(filename)
        except pd.io.common.EmptyDataError:
            pass


if __name__ == "__main__":
    main()
