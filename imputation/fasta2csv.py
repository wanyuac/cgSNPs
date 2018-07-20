#!/usr/bin/env python

"""
Generating a CSV-formatted SNP table from a multi-FASTA file.

Python 2 and 3 compatible
Copyright 2018 Yu Wan <wanyuac@gmail.com>
Licensed under the Apache License, Version 2.0
First and the latest edition: 19/7/2018
"""

from __future__ import print_function
import os
import sys
import textwrap
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import ArgumentParser, RawDescriptionHelpFormatter

prog_descr = """
Convert a multi-sequence alignment (in FASTA format) into a comma-delimited SNP table that is
compatible with the script parseSNPtable.py in the RedDog package. It is useful in processing
the SNPs after imputation.

Example: python fasta2csv -c snp.coords -i mapping/alignment_imputed.mfasta -o snps/snp_table.csv

Given an um-imputed SNP table from parseSNPtable.py, we can use the command
'cut -d, -f1 < snpTable.csv | tail -n +2 > snp.coords' to get the coordinate list while skipping
the first line of the column name 'Pos'.
  
Warning: this script requires sufficient memory to run when the FASTA file has a tremendous size.
"""


def parse_arguments():
    parser = ArgumentParser(prog = "fasta2csv", 
    formatter_class = RawDescriptionHelpFormatter,
    description = textwrap.dedent(prog_descr))
    parser.add_argument("-c", "--coords", dest = "coords", type = str, required = True,
                        help = "A single-column text file showing a list of SNP coordinates")
    parser.add_argument("-i", "--in", dest = "in_fasta", type = str, required = True,
                        help = "Input multi-FASTA file containing the sequence alignment")
    parser.add_argument("-k", "--key", dest = "key", type = str, required = False, default = "",
                        help = "A key word for strain names to be ignored.")  # useful for imputed SNP tables from ClonalFrameML, which contain sequences of internal nodes.
    parser.add_argument("-s", "--strains", dest = "strains", type = str, required = False, default = "",
                        help = "A single-column text file showing a list of strain names to be included.")  # This option applies after the 'key' option.
    parser.add_argument("-o", "--out", dest = "out", type = str, required = False, default = "snpTable.csv",
                        help = "Output SNP table")
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    coords = get_list(args.coords)
    if args.strains != "":
        ingroup = get_list(args.strains)
    else:
        ingroup = []
    out_csv = open(args.out, "w")  # exit if the file path is not accessible so that we can avoid the time-consuming step of fasta2dict
    snp_table = fasta2dict(args.in_fasta, coords, args.key, ingroup)  # may be memory consuming
    write_csv(snp_table, coords, out_csv)
    out_csv.close()
    print("The SNP table has been generated successfully.")
    
    return


def fasta2dict(mfasta, coords, key, ingroup):
    print("Importing sequence alignments and converting them into a matrix.")
    snp_tab = OrderedDict()
    n_coords = len(coords)
    key_filter = key != ""
    ingroup_filter = ingroup != []
    with open(mfasta, "rU") as f:
        for entry in SeqIO.parse(f, "fasta"):
            strain = entry.id
            
            # apply two levels of filter for strain names
            if key_filter:  # This scenario is more often than the ingroup filter in practice.
                if key in strain:  # For ClonalFrameML outputs, the key is "NODE_".
                    continue
            if ingroup_filter:  # apply a filter for sequence IDs
                if strain not in ingroup:
                    continue
            
            # store sequences in a dictionary
            seq = entry.seq  # a character string
            seq_bp = len(seq)
            check_lengths(n_coords, seq_bp, strain)  # Both lengths must equal.
            snp_tab[strain] = seq
            
    return snp_tab


def write_csv(tab, coords, out):
    # The argument "out" is a file handle.
    print("Writing the SNP table into a CSV file.")
    
    # print the header line
    strains = list(tab.keys())
    pos_num = len(coords)  # number of SNP positions
    out.write(",".join(["Pos"] + strains) + "\n")  # Pos,Strain1,Strain2,...

    # print SNPs at each position across strains
    for i in range(0, pos_num):
        line = [coords[i]]  # initialise this line using the current position
        for s in strains:
            line += tab[s][i]
        out.write(",".join(line) + "\n")

    return


def get_list(c):
    with open(c, "rU") as f:
        lst = f.read().splitlines()
    
    return lst


def check_lengths(n_coords, n_seq, strain):
    if n_coords != n_seq:
        print("Error: the number of SNP positions does not equal the alignment length of the strain %s." % strain)
        raise SystemExit
    
    return


if __name__ == "__main__":
    main()
