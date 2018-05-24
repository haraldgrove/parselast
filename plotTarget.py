#! /usr/bin/env python

"""
Create a plot of LAST alignment
Input file must contain only one target sequence.
Hits are displayed as boxes with colour based on name of query sequence
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.patches as patches


def roundup(x, n):
    N = pow(10, n)
    return x if x % N == 0 else x + N - x % N


def plot_coverage(filename):
    """
    Plots the location of hits
    """
    print("Reading {}".format(filename))
    df = pd.read_table(
        filename, comment="#", header=0, delim_whitespace=True, error_bad_lines=False
    )
    df["end1"] = df["start1"] + df["alnSize1"]
    df.sort_values(by=("start1"), inplace=True)
    minx, maxx = df["start1"].min(), df["end1"].max()
    unique = df["name2"].value_counts()
    fig, axes = plt.subplots(1, 1, figsize=(20, 20))
    axes.set_xlim([minx, maxx])
    axes.set_ylim([0, len(df)])
    norm = colors.Normalize(0, len(unique))
    rank = {}
    y = 0
    h = 1
    for (index, row) in df.iterrows():
        if row["qual"] != 1:
            # continue
            pass
        if row["name2"] not in rank:
            rank[row["name2"]] = len(rank)
        c = rank[row["name2"]]
        start = row["start1"]
        end = row["start1"] + row["alnSize1"]
        color = cm.viridis(norm(c))
        if row["qual"] == 1:
            axes.add_patch(
                patches.Rectangle((start, y), end - start, h, facecolor=color)
            )
        else:
            axes.add_patch(
                patches.Rectangle(
                    (start, y), end - start, h, facecolor=color, hatch="*"
                )
            )
        axes.text(
            start * 1.0001,
            y + 0.1,
            row["name2"],
            bbox={"facecolor": "white", "alpha": 1, "pad": 1},
        )
        y += h
    scale_s = "{:.0f}".format(0.1 * (maxx - minx))
    scale = roundup(int(scale_s), len(scale_s) - 1)
    axes.add_patch(
        patches.Rectangle((maxx - 2 * scale, 0 + 0.25), scale, 0.5, facecolor="black")
    )
    axes.text(
        maxx - 1.9 * scale, 0 + 0.3, "{} bp".format(scale), color="white", size=18
    )
    figurefile = filename.rsplit(".", 1)[0] + ".png"
    fig.savefig(figurefile, dpi=90, bbox_inches="tight")


def main():
    for filename in sys.argv[1:]:
        try:
            plot_coverage(filename)
        except pd.io.common.EmptyDataError:
            pass


if __name__ == "__main__":
    main()
