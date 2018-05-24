#!/usr/bin/env python3
"""
Module Docstring
"""

__author__ = "Harald Grove"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import time
import sys


def read_last(infile):
    db = {}
    with open(infile, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if name2 not in db:
                db[name2] = {"name1": name1}


def single_hit_qc(lastfile):
    """ Removes single lines from alignment based on:
        expected score vs. actual score
        adjusted for alignment size
    """
    db = {}
    with open(lastfile, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if int(alnSize2) == int(seqSize2):
                db[name2] = score
            if int(alnSize1) == int(seqSize1):
                db[name1] = score
    outfile1 = "{}.qc.txt".format(lastfile.rsplit(".", 1)[0])
    outfile2 = "{}.noise.txt".format(lastfile.rsplit(".", 1)[0])
    with open(lastfile, "r") as fin, open(outfile1, "w") as fout1, open(
        outfile2, "w"
    ) as fout2:
        for line in fin:
            if line.startswith("#"):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if int(alnSize2) < 1000:
                ratio = 6
            else:
                ratio = 5
            if int(score) < int(alnSize2) * ratio:
                fout2.write(line)
            else:
                if int(alnSize2) < int(seqSize2) and name2 in db:
                    fout2.write(line)
                elif int(alnSize1) < int(seqSize1) and name1 in db:
                    fout2.write(line)
                else:
                    fout1.write(line)


def read_graph(lastfile):
    """ Records Last entries in a dict:
        Key: sequence name (both name1 and name2)
        Value: Score and name of the highest scoring match
    """
    db = {}
    with open(lastfile, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if name2 not in db:
                db[name2] = {"score": int(score), "name1": name1, "links": []}
            else:
                if int(score) > db[name2]["score"]:
                    db[name2] = {"score": int(score), "name1": name1, "links": []}
            if name1 not in db:
                db[name1] = {"score": int(score), "name2": name2, "links": []}
            else:
                if int(score) > db[name1]["score"]:
                    db[name1] = {"score": int(score), "name2": name2, "links": []}
        for name, entry in db.items():
            if "name1" in entry:
                name1 = entry["name1"]
                db[name1]["links"].append(name)
            elif "name2" in entry:
                name2 = entry["name2"]
                db[name2]["links"].append(name)
    return db


def write_graph(db, outfile):
    """ Uses dict from 'read_graph' to output the best match from name1 to name2,
        and name2s best match from name1
    """
    with open(outfile, "w") as fout:
        for name, entry in db.items():
            if "name1" in entry:
                continue
            name1, score, name2 = name, entry["score"], entry["name2"]
            try:
                score2, oname1 = db[name2]["score"], db[name2]["name1"]
            except KeyError:
                print(name, name2)
                print(db[name2])
            fout.write(
                "{}\t{}\t{}\t{}\t{}\n".format(name1, score, name2, score2, oname1)
            )


def unique_match(lastfile):
    """ Identifies unique matches, i.e. sequences that only match to each other """
    db = {}
    with open(lastfile, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if name2 not in db:
                db[name2] = [name1]
            elif name1 not in db[name2]:
                db[name2].append(name1)
            if name1 not in db:
                db[name1] = [name2]
            elif name2 not in db[name1]:
                db[name1].append(name2)
    outfile1 = "{}.unique.txt".format(lastfile.rsplit(".", 1)[0])
    outfile2 = "{}.multi.txt".format(lastfile.rsplit(".", 1)[0])
    with open(outfile1, "w") as fout1, open(lastfile, "r") as fin, open(
        outfile2, "w"
    ) as fout2:
        for line in fin:
            if line.startswith("#"):
                fout1.write(line)
                fout2.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if len(db[name1]) == 1 and len(db[name2]) == 1 and db[name1][0] == name2:
                fout1.write(line)
            else:
                fout2.write(line)


def otm_match(lastfile):
    """ Identifies one-to-many matches """
    db = {}
    with open(lastfile, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if name2 not in db:
                db[name2] = [name1]
            elif name1 not in db[name2]:
                db[name2].append(name1)
            if name1 not in db:
                db[name1] = [name2]
            elif name2 not in db[name1]:
                db[name1].append(name2)
    outfile1a = "{}.otm1.txt".format(lastfile.rsplit(".", 1)[0])
    outfile1b = "{}.otm2.txt".format(lastfile.rsplit(".", 1)[0])
    outfile2 = "{}.mtm.txt".format(lastfile.rsplit(".", 1)[0])
    with open(lastfile, "r") as fin, open(outfile1a, "w") as fout1a, open(
        outfile1b, "w"
    ) as fout1b, open(outfile2, "w") as fout2:
        for line in fin:
            if line.startswith("#"):
                fout1a.write(line)
                fout1b.write(line)
                fout2.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if len(db[name1]) > 1 and len(db[name2]) > 1:
                fout2.write(line)
            elif len(db[name1]) > 1:
                for name in db[name1]:
                    if len(db[name]) > 1:
                        fout2.write(line)
                        break
                else:
                    fout1a.write(line)
            elif len(db[name2]) > 1:
                for name in db[name2]:
                    if len(db[name]) > 1:
                        fout2.write(line)
                        break
                else:
                    fout1b.write(line)


def repeat_ctgs(lastfile):
    """ Removes segments with large amount of hits """
    db = {}
    hits = {}
    with open(lastfile, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            # Prepare start&end position for both sequences
            start1 = int(start1)
            end1 = start1 + int(alnSize1)
            if strand2 == "+":
                start2x = int(start2)
                end2x = start2x + int(alnSize2)
            else:
                start2x = int(seqSize2) - (int(start2) + int(alnSize2))
                end2x = int(seqSize2) - (int(start2))
            # Add segment to list of already seen segments, update count of overlap
            if name2 not in db:
                db[name2] = [[start2x, end2x, 0]]
            else:
                count = 0
                ind = 0
                while ind < len(db[name2]):
                    s, e = db[name2][ind][0:2]
                    if (s + 100) < (end2x - 100) and (start2x + 100) < (e - 100):
                        count += 1
                        db[name2][ind][2] += 1
                    ind += 1
                db[name2].append([start2x, end2x, count])
            if name1 not in db:
                db[name1] = [[start1, end1, 0]]
            else:
                count = 0
                ind = 0
                while ind < len(db[name1]):
                    s, e = db[name1][ind][0:2]
                    if (s + 100) < (end1 - 100) and (start1 + 100) < (e - 100):
                        count += 1
                        db[name1][ind][2] += 1
                    ind += 1
                db[name1].append([start1, end1, count])
    outfile1 = "{}.normal.txt".format(lastfile.rsplit(".", 1)[0])
    outfile2 = "{}.repeat.txt".format(lastfile.rsplit(".", 1)[0])
    with open(lastfile, "r") as fin, open(outfile1, "w") as fout1, open(
        outfile2, "w"
    ) as fout2:
        for line in fin:
            if line.startswith("#"):
                fout1.write(line)
                fout2.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            rep = False
            for segment in db[name1]:
                if segment[2] > 1:
                    rep = True
                    break
            for segment in db[name2]:
                if segment[2] > 1:
                    rep = True
                    break
            if not rep:
                fout1.write(line)
            else:
                fout2.write(line)


def travel(name, db, ctgs, fout):
    for child in db[name]:
        try:
            ctgs.pop(child)
        except KeyError:
            continue
        fout.write("{}\n".format(child))
        travel(child, db, ctgs, fout)


def gather_network(lastfile):
    """ Identifies all connected sequences """
    db = {}
    hits = {}
    ctg_list = {}
    with open(lastfile, "r") as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if name2 not in db:
                db[name2] = [name1]
                hits[name2] = [line]
                ctg_list[name2] = 1
            elif name1 not in db[name2]:
                db[name2].append(name1)
                hits[name2].append(line)
            if name1 not in db:
                db[name1] = [name2]
                ctg_list[name1] = 1
            elif name2 not in db[name1]:
                db[name1].append(name2)
    outfile = "{}.groups.txt".format(lastfile.rsplit(".", 1)[0])
    gr = 0
    with open(outfile, "w") as fout:
        while len(ctg_list) > 0:
            name, value = ctg_list.popitem()
            fout.write("# group{}\n".format(gr))
            fout.write("{}\n".format(name))
            travel(name, db, ctg_list, fout)
            gr += 1


def filter_bed(lastfile, bedfile):
    """
    Selects lines from lastfile that overlaps the regions in bedfile
    :param lastfile:
    :param bedfile:
    :return:
    """
    db = {}
    with open(bedfile, "r") as fin:
        for line in fin:
            name, start, stop = line.strip().split()
            if name not in db:
                db[name] = []
            db[name].append([int(start), int(stop)])
    outfile1 = "{}.bed.txt".format(lastfile.rsplit(".", 1)[0])
    with open(lastfile, "r") as fin, open(outfile1, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if name1 not in db:
                continue
            # Prepare start&end position for both sequences
            start1 = int(start1)
            end1 = start1 + int(alnSize1)
            if strand2 == "+":
                start2x = int(start2)
                end2x = start2x + int(alnSize2)
            else:
                start2x = int(seqSize2) - (int(start2) + int(alnSize2))
                end2x = int(seqSize2) - (int(start2))
            for start, stop in db[name1]:
                if start1 < stop and start < end1:
                    fout.write(line)
                    break


def filter_gap_bed(lastfile, bedfile):
    """
        Selects lines from lastfile where a gap overlaps the regions in bedfile
        :param lastfile:
        :param bedfile:
        :return:
        """
    db = {}
    with open(bedfile, "r") as fin:
        for line in fin:
            name, start, stop, namex = line.strip().split()
            if name not in db:
                db[name] = []
            db[name].append([int(start), int(stop), namex])
    outfile1 = "{}.feature.txt".format(lastfile.rsplit(".", 1)[0])
    outfile2 = "target.bed"
    outfile3 = "query.bed"
    with open(lastfile, "r") as fin, open(outfile1, "w") as fout, open(
        outfile2, "w"
    ) as fout2, open(outfile3, "w") as fout3:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, name2, start2, alnSize2, strand2, seqSize2, blocks, *e = (
                line.strip().split()
            )
            if name1 not in db:
                continue
            # Prepare start&end position for both sequences
            start1 = int(start1)
            end1 = start1 + int(alnSize1)
            if strand2 == "+":
                start2x = int(start2)
                end2x = start2x + int(alnSize2)
            else:
                start2x = int(seqSize2) - (int(start2) + int(alnSize2))
                end2x = int(seqSize2) - (int(start2))
            gaps = []
            x1 = start1
            for b in blocks.split(","):
                if ":" in b:
                    s = x1
                    b1, b2 = b.split(":")
                    x1 += int(b1)
                    e = x1
                    gaps.append([s, e])
                else:
                    x1 += int(b)
            for start, stop, namex in db[name1]:
                for s1, e1 in gaps:
                    if s1 < stop and start < e1:
                        break
                else:
                    continue
                fout.write(line)
                fout2.write(
                    "{}\t{}\t{}\t{}_{}\n".format(
                        name1,
                        max(start - 100, 0),
                        min(stop + 100, int(seqSize1)),
                        name1,
                        namex,
                    )
                )
                shift = start - start1
                fout3.write(
                    "{}\t{}\t{}\t{}_{}\t0\t{}\n".format(
                        name2,
                        max(int(start2) + shift - 100, 0),
                        min(int(start2) + shift + (stop - start) + 100, int(seqSize2)),
                        name2,
                        namex,
                        strand2,
                    )
                )
                break


def main(args):
    """ Main entry point of the app """
    if args.option == "read_graph":
        db = read_graph(args.infile)
        write_graph(db, args.outfile)
    elif args.option == "single_hit_qc":
        single_hit_qc(args.infile)
    elif args.option == "unique_match":
        unique_match(args.infile)
    elif args.option == "otm_match":
        otm_match(args.infile)
    elif args.option == "repeat_ctgs":
        repeat_ctgs(args.infile)
    elif args.option == "gather_network":
        gather_network(args.infile)
    elif args.option == "filter_bed":
        filter_bed(args.infile, args.bedfile)
    elif args.option == "filter_gap_bed":
        filter_gap_bed(args.infile, args.bedfile)
    else:
        sys.stderr.write("Unknown operation: [{}]\n".format(args.option))
    if args.log:
        with open("README.txt", "a") as fout:
            fout.write("[{}]\t[{}]\n".format(time.asctime(), " ".join(sys.argv)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("option", help="Action to take")
    parser.add_argument("infile", help="Input file")

    # Optional argument flag which defaults to False
    parser.add_argument(
        "-l",
        "--log",
        action="store_true",
        default=False,
        help="Save command to 'README.txt'",
    )

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-o", "--outfile", help="Output file")
    parser.add_argument("-b", "--bedfile", help="bed-file")
    parser.add_argument("-n", "--name", action="store", dest="name")

    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v", "--verbose", action="count", default=0, help="Verbosity (-v, -vv, etc)"
    )

    # Specify output of '--version'
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__),
    )

    args = parser.parse_args()
    main(args)
