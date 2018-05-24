#! /usr/bin/env python

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.patches as patches


def locate_splits(target, query, gap):
    """
        Assigns chromosome location to each query contig
    """
    filename = "last{}_T{}_Q{}.E0.QC.txt".format(gap[0:3], target, query)
    names = [
        "score",
        "name1",
        "start1",
        "alnSize1",
        "strand1",
        "seqSize1",
        "name2",
        "start2",
        "alnSize2",
        "strand2",
        "seqSize2",
        "blocks",
        "EG",
        "E",
        "qual",
    ]
    df = pd.read_csv(filename, comment="#", header=0, delim_whitespace=True)
    fig, axes = plt.subplots(1, 1, figsize=(20, 20))
    split_contig = {}
    for (index, row) in df.iterrows():
        if row["qual"] != 1:
            continue
        # Only consider the Nuclear chromosomes
        try:
            rank = int(row["name1"][3:])
        except ValueError:
            continue
        if row["seqSize2"] < 1000:
            continue
        if row["name2"] not in split_contig:
            split_contig[row["name2"]] = {}
        if row["name1"] not in split_contig[row["name2"]]:
            split_contig[row["name2"]][row["name1"]] = 0
        split_contig[row["name2"]][row["name1"]] += row["alnSize2"]
    with open(
        "last{}_T{}_Q{}_splits.E0.QC.txt".format(gap[0:3], target, query), "w"
    ) as fout:
        for name2, hits in split_contig.items():
            for name1, length in hits.items():
                fout.write("{}\t{}\t{}\t{}\n".format(name2, name1, length, len(hits)))


def main(t, q, g):
    locate_splits(t, q, g)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        query = ["7002", "7011", "8001", "9011R", "9018R", "8006", "9023R"]
        target = "orsadb"
    else:
        query = [sys.argv[2]]
        target = sys.argv[1]
    gap = "gapped"
    for q in query:
        main(target, q, gap)
