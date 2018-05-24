#!/usr/bin/env python

import sys

chrom = sys.argv[1]
start, stop = int(sys.argv[2]), int(sys.argv[3])
for line in sys.stdin:
    if line.startswith("#"):
        continue
    line_l = line.strip().split()
    try:
        ch, st, en = line_l[1], int(line_l[2]), int(line_l[2]) + int(line_l[3])
    except ValueError:
        sys.stdout.write(line)
        continue
    if ch != chrom:
        continue
    if start < en and stop > st:
        sys.stdout.write(line)
