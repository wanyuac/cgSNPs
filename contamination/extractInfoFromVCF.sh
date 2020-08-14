#!/bin/bash
# Extract columns from a list of VCF files.
# Command line:
#   bash extractInfoFromVCF -i [a text file of sample names] -s [source directory] -o [output directory] \
#   -r [a common root (between the sample id and file name extension) shared by names of VCF files] \
#   -f [comma-delimited fields of INFO to be extracted] -d [directory of vcftools]
# Directory names must not be followed by forward lashes.
# The INFO domain must not contain POS. The position information is however always be extracted to get the cut command execute correctly.
# Example:
#   bash extractInfoFromVCF -i strains_all.txt -s vcf -o info -r '_study1' -f 'DP,DP4'
# Copyright (C) 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the GNU General Public License (GPL) version 3
# Earliest editions: 7, 20, 22/4/2017, 29/11/2017; the latest edition: 13 Aug 2020

vcftools_dir='.'  # Default path to vcftools is the current working directory, which is often not the case in the real world.

# Read five arguments
while [[ $# -gt 1 ]]; do  # loops as long as the number of arguments >= 2
    key="$1"
    case $key in
        -i|--input)
        id_file="$2"
        shift 2 # remove the first two arguments from the vector
        ;;
        -s|--source)
        source_dir="$2"
        shift 2
        ;;
        -o|--output)
        output_dir="$2"
        shift 2
        ;;
        -r|--root)
        filename_root="$2"
        shift 2
        ;;
        -f|--fields)
        IFS=',' read -r -a fields <<< "$2"  # convert a comma-delimited string into an array
        shift 2
        ;;
        -d|--vcftoolsDir)  # when it has been supplied in the command line, otherwise, use the default
        vcftools_dir="$2"
        shift 2
        ;;
        *)
        shift 2  # deliberately left blank for undefined arguments
        ;;
    esac  # end of the case statement
done

# Set up the output directory
if [ ! -d $output_dir ]; then
    mkdir $output_dir
fi

# Arrange fields to be extracted from VCFs
info=''
for v in "${fields[@]}"; do
    info="${info}--get-INFO ${v} "  # concatenate the "--get-INFO" arguments
done

# Four reserved columns in the output of vcftools: CHROM, POS, REF, ALT. They must not appear in the argument -f.
#n=$[ ${#fields[@]} + 5 ]  # POS, REF, ALT, FIELD_1 (the sixth column), FIELD_2, ... So it is 6 + n - 1 = 5 + n.

# Parse VCFs
while read -r id; do  # read lines in a file
	#${vcftools_dir}/vcftools --vcf ${source_dir}/${id}${filename_root}.vcf --out ${id} ${info} -c | cut -f2-4,6-${n} > ${output_dir}/${id}__info.tsv
    ${vcftools_dir}/vcftools --vcf ${source_dir}/${id}${filename_root}.vcf ${info}--stdout > ${output_dir}/${id}__info.tsv
done < "$id_file"
