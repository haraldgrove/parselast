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
    with open(infile, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                continue
            score, name1, start1, alnSize1, strand1, seqSize1, \
            name2, start2, alnSize2, strand2, seqSize2, blocks, *e =\
            line.strip().split()
            if name2 not in db:
                db[name2] = {'name1':name1}

def read_graph(lastfile):
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
            fout.write('{}\t{}\t{}\t{}\t{}\n'.format(name1,score,name2,score,oname1))


def main(args):
    """ Main entry point of the app """
    db = read_graph(args.infile)
    write_graph(db, args.outfile)
    if args.log:
        with open('README.txt', 'a') as fout:
            fout.write('[{}]\t[{}]\n'.format(time.asctime(), ' '.join(sys.argv)))


if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("infile", help="Input file")
    parser.add_argument("outfile", help="Output file")

    # Optional argument flag which defaults to False
    parser.add_argument('-l', '--log', action="store_true", default=False, help="Save command to 'README.txt'")

    # Optional argument which requires a parameter (eg. -d test)
    parser.add_argument("-n", "--name", action="store", dest="name")

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
