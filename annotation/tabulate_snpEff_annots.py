#!/usr/bin/env python

"""
Parse VCF/BCF-formatted output of snpEff and generates a tab-delimited file. This script only extracts
the ANN sub-domain from the INFO domain of the input file.

Command line: python tabulate_snpEff_annots.py -i annotated_variants.vcf -o parsed_annots.tsv

Dependencies: Python v3, package pysam (pysam.readthedocs.io).

Copyright 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 15 Apr 2020; latest modification: 29 Apr 2020
"""

from argparse import ArgumentParser
from pysam import VariantFile

"""
The ANN sub-domain in the INFO domain consists of 16 fields.
RANK records indices of annotations at each position as snpEff sorts the annotations based on an algorithm.
See snpeff.sourceforge.net/SnpEff_manual.html#input for details.
"""
HEADER = ['CHROM', 'POS', 'REF', 'ALT', 'RANK', 'EFFECT', 'IMPACT', 'GENE', 'LOCUS_TAG', 'FEATURE_TYPE',\
          'FEATURE_ID', 'TRANSCRIPT_BIOTYPE', 'RANK_TOTAL', 'NUCLEOTIDE_VAR', 'PROTEIN_VAR', 'CDNA_POS', 'CDS_POS',\
          'PROTEIN_POS', 'DISTANCE_TO_FEATURE', 'NOTE']


def parse_arguments():
    parser = ArgumentParser(description= 'Concatenate MLST alleles according to ST profiles')
    parser.add_argument('-i', type = str, required = True, help = 'Path to the input VCF/BCF')
    parser.add_argument('-o', type = str, required = True, help = 'Path to the output TSV')

    return parser.parse_args()


def main():
    args = parse_arguments()
    
    # Read the VCF/BCF file
    try:
        vcf = VariantFile(args.i)  # Automatic detection of file format
    except:
        print('The input file is not accessible.')
        
    # Parse the ANN sub-domain and write lines into the output file
    tsv = open(args.o, 'w')
    tsv.write('\t'.join(HEADER) + '\n')  # Write the header
    for rec in vcf.fetch():
        parseVCFRec(rec, tsv)
    tsv.close()
    
    return


def parseVCFRec(rec, tsv_handle):
    """
    Parse an object of the class pysam.libcbcf.VariantRecord.
    Use dir(rec) to see attributes of the record.
    """
    fields = [rec.chrom, rec.pos, rec.ref]
    annots = rec.info['ANN']  # A tuple object. Use len(annots) to show the number of its elements.
    i = 1  # The rank
    for f in annots:
        ann_fields = f.split('|')
        ann_fields = ['NA' if x == '' else x for x in ann_fields]  # Replace empty strings with 'NA'
        ann_fields = [rec.chrom, str(rec.pos), rec.ref, ann_fields[0], str(i)] + ann_fields[1:]  # Insert the rank into this list
        tsv_handle.write('\t'.join(ann_fields) + '\n')  # Write a new line into the output file
        i += 1  # Move to the annotation of the next rank
    
    return


if __name__ == '__main__':
    main()
