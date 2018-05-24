#! /usr/bin/env python

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.patches as patches


def prepare_last(target, query, gap):
    """
    Filter out hits with E-value above 0
    """
    header_l = [
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
    ]
    filename = "last{}_T{}_Q{}.txt".format(gap, target, query)
    df = pd.read_csv(
        filename, comment="#", header=None, names=header_l, delim_whitespace=True
    )
    # Remove hits with E-value different from 0 and write new file
    df1 = df[df["E"] == "E=0"]
    df1.to_csv(filename.rsplit(".", 1)[0] + ".E0.txt", sep="\t", index=None)


def filter_overlap(target, query, gap):
    """
        Indicate hits that are covered by other, higher scoring hits.
        Adds an extra column, "qual". 0 = removed, 1 = retained
    """
    filename = "last{}_T{}_Q{}.E0.txt".format(gap, target, query)
    df = pd.read_csv(filename, comment="#", header=0, delim_whitespace=True)
    split_contig = {}
    df["qual"] = np.repeat(1, len(df))
    for (index, row) in df.iterrows():
        # Only consider the Nuclear chromosomes
        try:
            rank = int(row["name1"][3:])
        except ValueError:
            row["qual"] = -1
            continue
        df1 = df[df["name2"] == row["name2"]]
        zstart, zend = row["start2"], row["start2"] + row["alnSize2"]
        if row["strand2"] == "-":
            zstart, zend = row["seqSize2"] - zend, row["seqSize2"] - zstart
        # if current line is lower score and either completely covered, or partly covered and on different chromosome
        for i1, r1 in df1.iterrows():
            zs, ze = r1["start2"], r1["start2"] + r1["alnSize2"]
            if r1["strand2"] == "-":
                zs, ze = r1["seqSize2"] - ze, r1["seqSize2"] - zs
            if row["score"] >= r1["score"]:
                continue
            if row["name1"] != r1["name1"]:
                if zs < zend and ze > zstart:
                    df.loc[index, "qual"] = 0
            else:
                if zstart >= zs and zend <= ze:
                    df.loc[index, "qual"] = 0
    df.to_csv(filename.rsplit(".", 1)[0] + ".QC.txt", sep="\t", index=None)


def main(t, q, g):
    prepare_last(t, q, g)
    filter_overlap(t, q, g)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        query = [
            "7002",
            "7003",
            "7010",
            "7011",
            "7012",
            "7016",
            "8001",
            "8002",
            "8003",
            "8005R",
            "8006",
            "9006R",
            "9007R",
            "9011R",
            "9012",
            "9018R",
            "9016R",
            "9023R",
        ]
        target = "orsadb"
    elif len(sys.argv) < 3:
        query = [
            "7002",
            "7003",
            "7010",
            "7011",
            "7012",
            "7016",
            "8001",
            "8002",
            "8003",
            "8005R",
            "8006",
            "9006R",
            "9007R",
            "9011R",
            "9012",
            "9018R",
            "9016R",
            "9023R",
        ]
        target = sys.argv[1]
    else:
        query = [sys.argv[2]]
        target = sys.argv[1]
    gap = "gap"
    for q in query:
        try:
            main(target, q, gap)
        except:
            pass
