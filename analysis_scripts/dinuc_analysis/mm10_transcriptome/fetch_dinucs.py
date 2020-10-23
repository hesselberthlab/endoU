#! /usr/bin/env python3

# module load python/3.8.5

from pyfaidx import Fasta, FetchError
import gzip
import pdb

def main(args):

    fasta = Fasta(args.fasta)
    tab = gzip.open(args.table)

    for line in tab:
        name, pos, count, norm_count = [f.strip() for f in line.decode().split("\t")]

        # depth files start at 1, fasta at 0
        pos = int(pos) - 1

        start = pos + args.up
        end = pos + args.down

        if start >= end: continue

        try:
            dinuc = fasta[name][start:end].seq
        except FetchError:
            continue  

        # write the primary id
        id = name.split("|")[0]

        print(id, pos, count, norm_count,  dinuc, sep = "\t")

if __name__ == "__main__":
   
    from argparse import ArgumentParser

    parser = ArgumentParser("get dinucleotides")

    parser.add_argument("fasta", help="fasta file")
    parser.add_argument("table", help="table")

    parser.add_argument("--up", type=int, help="upstream offset (can be negative)")
    parser.add_argument("--down", type=int, help="downstream offset")

    args = parser.parse_args()

    main(args)
