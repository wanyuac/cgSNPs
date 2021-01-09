#!/usr/bin/env python

"""
Generating a CSV-formatted SNP table from a multi-FASTA file. The output CSV file can be processed
by parseSNPtable.py of package RedDog and function importCoreGenomeSNPs of R package GeneMates.

Note that the easiest way for converting a TSV file from Snippy to a GeneMates-compatible CSV file
is using the command: sed 's/\t/,/g' core.tab | cut -d',' -f2- > core.csv rather than this script.

Python 2 and 3 compatible
Copyright 2018 Yu Wan <wanyuac@gmail.com>
Licensed under the Apache License, Version 2.0
First edition: 19-20 July 2018; the latest edition: 9 Jan 2021
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
compatible with the script parseSNPtable.py in the RedDog package as well as function
importCoreGenomeSNPs of R package GeneMates. It is useful in processing SNPs after imputation
and also in visualising SNPs as a integer matrix. This script assumes equal lengths of all sequences
and the script terminates otherwise.

Example: python fasta2csv.py -i mapping/alignment_imputed.mfasta -o snps/snp_table.csv -c snp.coords

Given an um-imputed SNP table from parseSNPtable.py, we can use the command
'cut -d, -f1 < snpTable.csv | tail -n +2 > snp.coords' to get the coordinate list while skipping
the first line of the column name 'Pos'.
  
Warning: this script requires sufficient memory to run when the FASTA file has a tremendous size.
"""


def parse_arguments():
    parser = ArgumentParser(prog = "fasta2csv", 
    formatter_class = RawDescriptionHelpFormatter,
    description = textwrap.dedent(prog_descr))
    parser.add_argument("-i", "--in", dest = "in_fasta", type = str, required = True,
                        help = "Input multi-FASTA file containing the sequence alignment")
    parser.add_argument("-c", "--coords", dest = "coords", type = str, required = False, default = "",
                        help = "An optional single-column text file showing a list of SNP coordinates")
    parser.add_argument("-k", "--key", dest = "key", type = str, required = False, default = "",
                        help = "An optional key word for strain names to be ignored.")  # useful for imputed SNP tables from ClonalFrameML, which contain sequences of internal nodes.
    parser.add_argument("-s", "--strains", dest = "strains", type = str, required = False, default = "",
                        help = "An optional single-column text file showing a list of strain names to be included.")  # This option applies after the 'key' option.
    parser.add_argument("-o", "--out", dest = "out", type = str, required = False, default = "snpTable.csv",
                        help = "Output SNP table")
    
    return parser.parse_args()


def main():
    args = parse_arguments()
    coords = get_list(args.coords)
    ingroup = get_list(args.strains)
    out_csv = open(args.out, "w")  # exit if the file path is not accessible so that we can avoid the time-consuming step of fasta2dict
    snp_table = fasta2dict(args.in_fasta, coords, args.key, ingroup)  # may be memory consuming
    write_csv(snp_table, coords, out_csv)
    out_csv.close()
    print("The SNP table has been generated successfully.")
    
    return


def fasta2dict(mfasta, coords, key, ingroup):
    print("Importing sequence alignments and converting them into a matrix.")
    #snp_tab = OrderedDict()
    snp_tab = dict()
    use_arbitrary_coords = (coords == None)
    seq_len_ref = len(coords) if not use_arbitrary_coords else 0  # Reference length of sequences 
    key_filter = (key != "")
    ingroup_filter = (ingroup != None)
    if (not key_filter) and (not ingroup_filter):
        print("No filter for sequences is applied.")
    
    with open(mfasta, "r") as f:
        for entry in SeqIO.parse(f, "fasta"):
            strain = entry.id
            
            # apply two levels of filter for strain names
            if key_filter:  # This scenario is more often than the ingroup filter in practice.
                if key in strain:  # For ClonalFrameML outputs, the key is "NODE_".
                    continue
            if ingroup_filter:  # apply a filter for sequence IDs
                if strain not in ingroup:
                    continue
            
            # store sequences in dictionary snp_tab
            seq = entry.seq  # a character string
            seq_bp = len(seq)
            if use_arbitrary_coords:
                if seq_len_ref > 0:
                    check_lengths(seq_len_ref, seq_bp, strain)  # Both lengths must be equal.
                else:
                    seq_len_ref = seq_bp  # Set seq_len_ref to the length of the first input sequence
            else:
                check_lengths(seq_len_ref, seq_bp, strain)
            snp_tab[strain] = seq
            
    return snp_tab


def write_csv(tab, coords, out):
    # The argument "out" is a file handle.
    print("Writing the SNP table into a CSV file.")
    
    # print the header line
    strains = list(tab.keys())
    out.write(",".join(["Pos"] + strains) + "\n")  # Pos,Strain1,Strain2,...
    
    # Create arbitary SNP coordinates when it is necessary
    if coords == None:
        coords = [str(i) for i in range(1, len(tab[strains[0]]) + 1)]
    pos_num = len(coords)  # number of SNP positions

    # print SNPs at each position across strains
    for i in range(0, pos_num):
        line = [coords[i]]  # initialise this line using the current position
        for s in strains:
            line += tab[s][i]  # Extend the list
        out.write(",".join(line) + "\n")

    return


def get_list(c):
    if c != "":
        if os.path.exists(c):
            with open(c, "r") as f:
                lst = f.read().splitlines()
        else: lst = None
    else:
        lst = None
    
    return lst


def check_lengths(n_coords, n_seq, strain):
    if n_coords != n_seq:
        print("Error: the number of SNP positions does not equal the alignment length of the strain %s." % strain)
        raise SystemExit
    
    return


if __name__ == "__main__":
    main()
