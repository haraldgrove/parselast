#!/usr/bin/env python
import sys


def convertmaf(maf):
    names = []
    db = {}
    gapdb = {}
    ref = True
    with open(maf, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            line_l = line.strip().split()
            if len(line_l) < 7:
                continue
            if line_l[0] == "a":
                continue
            code1, name1, start1, end1, strand1, size1, seq1 = line_l
            code2, name2, start2, end2, strand2, size2, seq2 = next(fin).strip().split()
            name2 = name2.split("_")[0]
            if name1 not in names:
                names.append(name1)
            if name2 not in names:
                names.append(name2)
            pos = int(start1)
            gap = 0
            gapseq = ""
            for ind, base1 in enumerate(seq1):
                if base1 != "-":
                    if len(gapseq) > 0:
                        gapdb[pos][name2] = gapseq
                    gap, gapseq = 0, ""
                    pos += 1
                    if pos not in db:
                        db[pos] = {}
                    db[pos][name1] = base1
                    if name2 in db[pos]:
                        if db[pos][name2] != seq2[ind]:
                            db[pos][name2] = "X"
                    else:
                        db[pos][name2] = seq2[ind]
                else:
                    if pos not in gapdb:
                        gapdb[pos] = {}
                    gapseq += seq2[ind]
    return db, gapdb, int(size1), names


def convertmaf2(maf):
    names = []
    db = {}
    gapdb = {}
    ref = True
    with open(maf, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            line_l = line.strip().split()
            if len(line_l) < 7:
                continue
            if line_l[0] == "a":
                continue
            code1, name1, start1, end1, strand1, size1, seq1 = line_l
            code2, name2, start2, end2, strand2, size2, seq2 = next(fin).strip().split()
            name2 = name2.split("_")[0]
            if name1 not in names:
                names.append(name1)
            if name2 not in names:
                names.append(name2)
            pos = int(start1)
            gap = 0
            gapseq = ""
            for ind, base1 in enumerate(seq1):
                if base1 != "-":
                    if len(gapseq) > 0:
                        gapdb[pos][name2] = gapseq
                    gap, gapseq = 0, ""
                    pos += 1
                    if pos not in db:
                        db[pos] = {}
                    db[pos][name1] = base1
                    if name2 in db[pos]:
                        if db[pos][name2] != seq2[ind]:
                            db[pos][name2] = "X"
                    else:
                        db[pos][name2] = seq2[ind]
                else:
                    if pos not in gapdb:
                        gapdb[pos] = {}
                    gapseq += seq2[ind]
    return db, gapdb, int(size1), names


def printmaf(db, gapdb, size, names):
    out = {}
    for name in names:
        out[name] = ""
    for pos in range(0, size):
        try:
            seq = db[pos]
        except KeyError:
            for name in names:
                out[name] += "-"
            continue
        for name in names:
            try:
                base = seq[name]
            except KeyError:
                base = "-"
            out[name] += base
        if pos in gapdb:
            gap = gapdb[pos]
            gaplen = 0
            # Find the longest inserted segment
            for name, gapseq in gap.items():
                gaplen = max(gaplen, len(gapseq))
            for name in out:
                try:
                    gapseq = gap[name]
                    out[name] += gapseq.ljust(gaplen, "-")
                except KeyError:
                    out[name] += "-" * gaplen
    for name, seq in out.items():
        sys.stdout.write(">{}\n{}\n".format(name, seq))


def main():
    maffile = sys.argv[1]
    # db, gapdb, size, names = convertmaf(maffile)
    db, gapdb, size, names = convertmaf2(maffile)
    printmaf(db, gapdb, size, names)


if __name__ == "__main__":
    main()
