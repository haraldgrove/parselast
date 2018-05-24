#!/usr/bin/env python

"""
Filter out hits where the contig from Query/Target is completely covered by a contig from Target/Query
Target is prioritized, and there have to be a substantial size difference
"""

import sys
import re

HEADER = (
    "score\tname1\tstart1\talnSize1\tstrand1\tseqSize1\tname2\tstart2\talnSize2\tstrand2\tseqSize2\tblocks\tEG\tE"
)


def split_lines(target, query):
    written = {}
    infile = sys.argv[1]
    out1 = "T{}.redundant.txt".format(target)
    out2 = "Q{}.redundant.txt".format(query)
    with open(infile, "r") as fin, open(out1, "w") as fout1, open(out2, "w") as fout2:
        fout1.write("Kept_Target\tRedundant_Query\n")
        fout2.write("Kept_Query\tRedundant_Target\n")
        for line in fin:
            if line.startswith("#"):
                continue
            line_l = line.strip().split()
            try:
                name1, start1, alnSize1, seqSize1 = (
                    line_l[1],
                    int(line_l[2]),
                    int(line_l[3]),
                    int(line_l[5]),
                )
                name2, start2, alnSize2, seqSize2 = (
                    line_l[6],
                    int(line_l[7]),
                    int(line_l[8]),
                    int(line_l[10]),
                )
            except ValueError:
                continue
            if (
                alnSize1 >= seqSize1
                and alnSize2 < seqSize2
                and seqSize1 < seqSize2 * 0.9
            ):
                fout2.write("{}\t{}\n".format(name2, name1))
            if (
                alnSize2 >= seqSize2
                and alnSize1 < seqSize1
                and seqSize2 < seqSize1 * 0.9
            ):
                fout1.write("{}\t{}\n".format(name1, name2))


def main():
    target = re.search("T([A-Za-z0-9]*)", sys.argv[1])
    query = re.search("Q([A-Za-z0-9]*)", sys.argv[1])
    split_lines(target.group(1), query.group(1))


if __name__ == "__main__":
    main()
