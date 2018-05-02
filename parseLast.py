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
import pandas as pd


class Last(object):

    def __init__(self, lastfile):
        self.lastfile = lastfile
        self.lines = []
        self.header = []

    def _calc_seqid(self, score, blocks):
        alnSize = 0
        gaps = 0
        gapSize = 0
        for b in blocks.split(','):
            if ':' in b:
                gap = max([int(a) for a in b.split(':')])
                gaps += gap
                gapSize += (21 + 9 * gap)
            else:
                alnSize += int(b)
        mismatch = ((alnSize * 6 - gapSize) - score ) / 24
        return 1 - (mismatch + gaps) / (gaps + alnSize)

    def read_last(self, idlim = 0, covlim = 0, lenlim = 0):
        """
        Reads each alignment from the last-alignment file.
        :param makedb:
        :param idlim: Minimum percent identity, in percent
        :param covlim: Minimum coverage, in percent
        :param lenlim: Minimum alignment length
        :return:
        """
        head = True
        with open(self.lastfile, 'r') as fin:
            for line in fin:
                if line.startswith('#'):
                    if head:
                        self.header.append(line)
                    continue
                head = False
                l = line.strip().split()
                if len(l) != 14:
                    continue
                score, start1, alnSize1, seqSize1, \
                start2, alnSize2, seqSize2 = [int(i) for i in [l[0], l[2], l[3], l[5], l[7], l[8], l[10]]]
                name1, strand1, \
                name2, strand2, blocks, *e = [l[1], l[4], l[6], l[9], l[11],l[12:]]
                end1 = start1 + alnSize1
                end2 = start2 + alnSize2
                if strand2 == '-':
                    start2p = seqSize2 - end2
                    end2p = seqSize2 - start2
                else:
                    start2p = start2
                    end2p = end2
                seqid = 100 * self._calc_seqid(score, blocks)
                if seqSize2 <= seqSize1:
                    seqcov = 100 * (alnSize2 / seqSize2)
                else:
                    seqcov = 100 * (alnSize1 / seqSize1)
                # Deciding to keep the alignment
                if seqid < idlim:
                    continue
                if alnSize2 < lenlim:
                    continue
                if seqcov < covlim:
                    continue
                self.lines.append([score, seqid, seqcov, name1, start1, alnSize1, end1, strand1, seqSize1,
                                   name2, start2, alnSize2, end2, strand2, seqSize2, blocks,
                                   start2p, end2p])


    def write_last(self, fout = sys.stdout):
        for line in self.header:
            l = line.strip().split()
            if len(l) > 1 and l[1] == 'score':
                l.insert(2,'idpct')
                l.insert(3,'covpct')
                l.insert(7, 'end1')
                l.insert(13,'end2')
                l.append('start2+')
                l.append('end2+')
            else:
                continue
            line = '{}\n'.format('\t'.join(l[1:]))
            fout.write(line)
        for line in self.lines:
            output = '\t'.join([str(b) for b in line])
            fout.write('{}\n'.format(output))

class Last2(object):

    def __init__(self, lastfile):
        self.lastfile = lastfile
        self.lines = []
        self.db = {}
        self.header = []
        self.name1list = []
        self.name2list = []

    def read_last(self):
        """
        columns: 0 - score, 1 - seqid, 2 - seqcov,
                 3 - name1, 4 - start1, 5 - alnSize1, 6 - end1, 7 - strand1, 8 - seqSize1,
                 9 - name2, 10 - start2, 11 - alnSize2, 12 - end2, 13 - strand2, 14 - seqSize2, 15 - blocks,
                 16 - start2p, 17 - end2p
        :return:
        """
        self.df = pd.read_table(self.lastfile, header = 0)

    def check_name2(self, name):
        df1 = self.df[['name2']==name]
        df1.sort_values(by='start2p')
        for ind in range(1,len(df1)-1):
            row0 = df1[ind-1]
            row1 = df1[ind]
            row2 = df1[ind+1]
            before2 = row1['start2p'] - row0['end2p']
            after2 = row2['start2p'] - row1['end2p']
            if row0['name1'] == row1['name1']: # This part assumes that 'name1' is not misassembled...
                if row1['strand2'] == '+':
                    before1 = row1['start1'] - row0['end1']
                else:
                    before1 = row0['start1'] - row1['end1']
            else:
                if row1['strand2'] == '+':
                    before1_1 = row1['start1']
                else:
                    before1_1 = (row1['seqSize1'] - row1['end1'])
                if row0['strand2'] == '+':
                    before1_2 = (row0['seqSize1'] - row0['end1'])
                else:
                    before1_2 = row0['start1']
                before1 = before1_1 + before1_2
                if row1['strand2'] == '+':
                    after1_1 = (row1['seqSize1'] - row1['end1'])
                else:
                    after1_1 = row1['start1']
                if row2['strand2'] == '+':
                    after1_2 = row2['start1']
                else:
                    after1_2 = (row2['seqSize1'] - row2['end1'])
                after1 = after1_1 + after1_2
            print(ind, before1, before2)
            print(ind, after1, after2)

def single_hit_qc(lastfile):
    """ Removes single lines from alignment based on:
        expected score vs. actual score
        adjusted for alignment size
    """
    db = {}
    with open(lastfile, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e =\
            line.strip().split()
            if int(alnSize2) == int(seqSize2):
                db[name2] = score
            if int(alnSize1) == int(seqSize1):
                db[name1] = score
    outfile1 = '{}.qc.txt'.format(lastfile.rsplit('.', 1)[0])
    outfile2 = '{}.noise.txt'.format(lastfile.rsplit('.', 1)[0])
    with open(lastfile, 'r') as fin, open(outfile1, 'w') as fout1, open(outfile2, 'w') as fout2:
        for line in fin:
            if line.startswith('#'):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e =\
            line.strip().split()
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
    with open(lastfile, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
            if name2 not in db:
                db[name2] = {'score':int(score), 'name1':name1, 'links':[]}
            else:
                if int(score) > db[name2]['score']:
                    db[name2] = {'score': int(score), 'name1': name1, 'links':[]}
            if name1 not in db:
                db[name1] = {'score':int(score), 'name2':name2, 'links':[]}
            else:
                if int(score) > db[name1]['score']:
                    db[name1] = {'score': int(score), 'name2': name2, 'links':[]}
        for name, entry in db.items():
            if 'name1' in entry:
                name1 = entry['name1']
                db[name1]['links'].append(name)
            elif 'name2' in entry:
                name2 = entry['name2']
                db[name2]['links'].append(name)
    return db

def write_graph(db,outfile):
    """ Uses dict from 'read_graph' to output the best match from name1 to name2,
        and name2s best match from name1
    """
    with open(outfile,'w') as fout:
        for name,entry in db.items():
            if 'name1' in entry:
                continue
            name1, score, name2 = name, entry['score'], entry['name2']
            try:
                score2, oname1 = db[name2]['score'], db[name2]['name1']
            except KeyError:
                print(name,name2)
                print(db[name2])
            fout.write('{}\t{}\t{}\t{}\t{}\n'.format(name1,score,name2,score2,oname1))

def unique_match(lastfile):
    """ Identifies unique matches, i.e. sequences that only match to each other """
    db = {}
    with open(lastfile, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
            if name2 not in db:
                db[name2] = [name1]
            elif name1 not in db[name2]:
                db[name2].append(name1)
            if name1 not in db:
                db[name1] = [name2]
            elif name2 not in db[name1]:
                db[name1].append(name2)
    outfile1 = '{}.unique.txt'.format(lastfile.rsplit('.',1)[0])
    outfile2 = '{}.multi.txt'.format(lastfile.rsplit('.',1)[0])
    with open(outfile1,'w') as fout1, open(lastfile, 'r') as fin, open(outfile2, 'w') as fout2:
        for line in fin:
            if line.startswith('#'):
                fout1.write(line)
                fout2.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
            if len(db[name1]) == 1 and len(db[name2]) == 1 and db[name1][0] == name2:
                fout1.write(line)
            else:
                fout2.write(line)

def otm_match(lastfile):
    """ Identifies one-to-many matches """
    db = {}
    with open(lastfile, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
            if name2 not in db:
                db[name2] = [name1]
            elif name1 not in db[name2]:
                db[name2].append(name1)
            if name1 not in db:
                db[name1] = [name2]
            elif name2 not in db[name1]:
                db[name1].append(name2)
    outfile1a = '{}.otm1.txt'.format(lastfile.rsplit('.',1)[0])
    outfile1b = '{}.otm2.txt'.format(lastfile.rsplit('.',1)[0])
    outfile2 = '{}.mtm.txt'.format(lastfile.rsplit('.',1)[0])
    with open(lastfile, 'r') as fin, open(outfile1a,'w') as fout1a, open(outfile1b,'w') as fout1b, open(outfile2, 'w') as fout2:
        for line in fin:
            if line.startswith('#'):
                fout1a.write(line)
                fout1b.write(line)
                fout2.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
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
    with open(lastfile, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
            # Prepare start&end position for both sequences
            start1 = int(start1)
            end1 = start1 + int(alnSize1)
            if strand2 == '+':
                start2x = int(start2)
                end2x = start2x + int(alnSize2)
            else:
                start2x = int(seqSize2) - (int(start2) + int(alnSize2))
                end2x = int(seqSize2) - (int(start2))
            # Add segment to list of already seen segments, update count of overlap
            if name2 not in db:
                db[name2] = [[start2x,end2x,0]]
            else:
                count = 0
                ind = 0
                while ind < len(db[name2]):
                    s,e = db[name2][ind][0:2]
                    if (s+100) < (end2x-100) and (start2x+100) < (e-100):
                        count += 1
                        db[name2][ind][2] += 1
                    ind += 1
                db[name2].append([start2x,end2x,count])
            if name1 not in db:
                db[name1] = [[start1,end1,0]]
            else:
                count = 0
                ind = 0
                while ind < len(db[name1]):
                    s,e = db[name1][ind][0:2]
                    if (s+100) < (end1-100) and (start1+100) < (e-100):
                        count += 1
                        db[name1][ind][2] += 1
                    ind += 1
                db[name1].append([start1,end1,count])
    outfile1 = '{}.normal.txt'.format(lastfile.rsplit('.',1)[0])
    outfile2 = '{}.repeat.txt'.format(lastfile.rsplit('.',1)[0])
    with open(lastfile, 'r') as fin, open(outfile1,'w') as fout1, open(outfile2, 'w') as fout2:
        for line in fin:
            if line.startswith('#'):
                fout1.write(line)
                fout2.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
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

def travel(name,db, ctgs, fout):
    for child in db[name]:
        try:
            ctgs.pop(child)
        except KeyError:
            continue
        fout.write('{}\n'.format(child))
        travel(child,db,ctgs,fout)

def gather_network(lastfile):
    """ Identifies all connected sequences """
    db = {}
    hits = {}
    ctg_list = {}
    with open(lastfile, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
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
    outfile = '{}.groups.txt'.format(lastfile.rsplit('.', 1)[0])
    gr = 0
    with open(outfile, 'w') as fout:
        while len(ctg_list) > 0:
            name, value = ctg_list.popitem()
            fout.write('# group{}\n'.format(gr))
            fout.write('{}\n'.format(name))
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
    with open(bedfile, 'r') as fin:
        for line in fin:
            name, start, stop = line.strip().split()
            if name not in db:
                db[name] = []
            db[name].append([int(start),int(stop)])
    outfile1 = '{}.bed.txt'.format(lastfile.rsplit('.', 1)[0])
    with open(lastfile, 'r') as fin, open(outfile1, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
            if name1 not in db:
                continue
            # Prepare start&end position for both sequences
            start1 = int(start1)
            end1 = start1 + int(alnSize1)
            if strand2 == '+':
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
    with open(bedfile, 'r') as fin:
        for line in fin:
            name, start, stop, namex = line.strip().split()
            if name not in db:
                db[name] = []
            db[name].append([int(start), int(stop), namex])
    outfile1 = '{}.feature.txt'.format(lastfile.rsplit('.', 1)[0])
    outfile2 = 'target.bed'
    outfile3 = 'query.bed'
    with open(lastfile, 'r') as fin, open(outfile1, 'w') as fout, \
            open(outfile2, 'w') as fout2, open(outfile3, 'w') as fout3:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e = \
                line.strip().split()
            if name1 not in db:
                continue
            # Prepare start&end position for both sequences
            start1 = int(start1)
            end1 = start1 + int(alnSize1)
            if strand2 == '+':
                start2x = int(start2)
                end2x = start2x + int(alnSize2)
            else:
                start2x = int(seqSize2) - (int(start2) + int(alnSize2))
                end2x = int(seqSize2) - (int(start2))
            gaps = []
            x1 = start1
            for b in blocks.split(','):
                if ':' in b:
                    s = x1
                    b1,b2 = b.split(':')
                    x1 += int(b1)
                    e = x1
                    gaps.append([s,e])
                else:
                    x1 += int(b)
            for start, stop, namex in db[name1]:
                for s1, e1 in gaps:
                    if s1 < stop and start < e1:
                        break
                else:
                    continue
                fout.write(line)
                fout2.write('{}\t{}\t{}\t{}_{}\n'.format(name1, max(start-100,0), min(stop+100,int(seqSize1)), name1,namex))
                shift = start - start1
                fout3.write('{}\t{}\t{}\t{}_{}\t0\t{}\n'.
                            format(name2, max(int(start2)+shift-100, 0),
                                   min(int(start2)+shift+(stop-start)+100, int(seqSize2)), name2, namex, strand2))
                break

def main(args):
    """ Main entry point of the app """
    last = Last(args.infile)
    idlim, covlim, lenlim = [int(a) for a in args.limits.split(',')]
    last.read_last(idlim=idlim, covlim=covlim, lenlim=lenlim)
    if args.outfile is not None:
        last.write_last(open(args.outfile, 'w'))
    else:
        last.write_last()
    if args.log:
        with open('README.txt', 'a') as fout:
            fout.write('[{}]\t[{}]\n'.format(time.asctime(), ' '.join(sys.argv)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    #parser.add_argument("option", help="Action to take")
    parser.add_argument("infile", help="Input file")

    # Optional argument flag which defaults to False
    parser.add_argument('-l', '--log', action="store_true", default=False, help="Save command to 'README.txt'")

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-o", "--outfile", help="Output file")
    parser.add_argument("-c", "--limits", default="0,0,0", help="Minimum idpct, covpct and length")
    #parser.add_argument("-b", "--bedfile", help='bed-file')
    #parser.add_argument("-n", "--name", action="store", dest="name")

    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        '-v',
        '--verbose',
        action='count',
        default=0,
        help="Verbosity (-v, -vv, etc)")

    # Specify output of '--version'
    parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s (version {version})'.format(version=__version__))

    args = parser.parse_args()
    main(args)
