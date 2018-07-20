#!/usr/bin/env python

"""
Python 2 and 3 compatible
Copyright 2018 Yu Wan <wanyuac@gmail.com>
Licensed under the Apache License, Version 2.0
First and the latest edition: 20/7/2018
"""

from __future__ import print_function
import os
import sys
import textwrap
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser, RawDescriptionHelpFormatter

prog_descr = """
Converting SNP imputation outputs of ClonalFrameML v1.11 into a multi-FASTA file so that
we can use fasta2csv.py to generate an imputed SNP table for subsequent analyses.

Example: python clonalFrameML2Fasta.py -p demo.ML_sequence.fasta -m demo.position_cross_reference.txt -o demo.mfasta
"""


def parse_arguments():
    parser = ArgumentParser(prog = "fasta2csv", 
    formatter_class = RawDescriptionHelpFormatter,
    description = textwrap.dedent(prog_descr))
    parser.add_argument("-p", "--pats", dest = "pats", type = str, required = True,
                        help = "A multi-FASTA file for SNP patterns (the ML_sequence file)")
    parser.add_argument("-m", "--mapping", dest = "hash", type = str, required = True,
                        help = "A comma-delimited file (the position_cross_reference file) mapping SNPs to their patterns")
    parser.add_argument("-o", "--out", dest = "out", type = str, required = False, default = "imputed.mfasta",
                        help = "The output multi-FASTA file")
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    pos = read_hash(args.hash)  # mapping SNPs to their positions in the pattern file
    seqs = restore_seqs(pos, args.pats, args.out)
    
    return


def read_hash(f):
    with open(f, "rU") as m:
        pos = m.readline().strip("\n")
    pos_str = pos.split(",")
    pos = []
    for i in range(0, len(pos_str)):
        pos += [int(pos_str[i]) - 1]  # convert to the character position in a string for Python
    
    return pos


def restore_seqs(pos, pat_file, out_file):
    p = open(pat_file, "rU")
    o = open(out_file, "w")
    patterns = SeqIO.parse(p, "fasta")
    seq_len = len(pos)
    for pat_seq in patterns:
        strain = pat_seq.id
        pats = pat_seq.seq
        s = ""
        for i in range(0, seq_len):  # go through each position in the original sequence
            s += pats[pos[i]]
        o.write(">" + strain + "\n")
        o.write(s + "\n")            
    p.close()
    o.close()
    
    return


if __name__ == "__main__":
    main()
