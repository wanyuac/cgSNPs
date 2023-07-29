#!/bin/bash
# This pipeline maps paired-end Illumina reads and Nanopore long reads against a reference genome and produces BAM files for variant calling.
# Dependency: samtools v1.17 and later versions, because earlier versions do not have option '-o'.
# Copyright (C) 2023 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Release: 11 Jan 2023; latest update: 29 July 2023

# Help information ###############
show_help() {
    echo "
    Download files of paired-end short reads from the NCBI SRA database.
    Arguments:
        -1: File of forward Illumina reads (R1)
        -2: File of reverse Illumina reads (R2)
        -l: File of long reads
        -r: File of the reference sequence
        -s: Output sample name
        -t: Number of threads
    Example command:
        ./hybrid_read_mapping.sh -1=sample1_1.fastq.gz -2=sample1_2.fastq.gz -l=sample1.fastq.gz -r=reference.fna -s=sample1 -t=4
    Dependencies: minimap2, samtools
    "
}

if [ -z "$1" ] || [ $1 = "-h" ]; then
    show_help
    exit
fi

# Read arguments ###############
for i in "$@"; do
    case $i in
        -1=*)
        r1="${i#*=}"
        ;;
        -2=*)
        r2="${i#*=}"
        ;;
        -l=*)
        l="${i#*=}"
        ;;
        -r=*)
        r="${i#*=}"
        ;;
        -s=*)
        s="${i#*=}"
        ;;
        -t=*)
        t="${i#*=}"
        ;;
        *)  # Do nothing otherwise.
        ;;
    esac
done

# Print arguments for subsequent STEPS ===============
echo 'Arguments for this pipeline:'
echo "  Illumina reads: ${r1}, ${r2}"
echo "  Nanopore reads: $l"
echo "  Reference genome: $r"
echo "  Output sample name: $s"
echo "  Number of cores: $t"
echo '  Read mapper: minimap2'

# Read mapping ###############
samtools faidx $r

echo 'Mapping Illumina reads against the reference genome'
minimap2 --MD -a -x sr -t $t $r $r1 $r2 > ${s}_illumina.sam

echo 'Mapping Nanopore reads against the reference genome'
minimap2 --MD -a -x map-ont -t $t $r $l > ${s}_nanopore.sam

#Creating BAM files ###############

# Create position-sorted BAM files ===============
echo 'Creating BAM files sorted by positions'
samtools view -S -b -u -@ $t ${s}_illumina.sam | samtools sort -@ $t -o ${s}_illumina_sorted.bam -
samtools view -S -b -u -@ $t ${s}_nanopore.sam | samtools sort -@ $t -o ${s}_nanopore_sorted.bam -
if [ -f "${s}_illumina_sorted.bam" ] && [ -f "${s}_nanopore_sorted.bam" ]; then
    samtools index -@ $t ${s}_illumina_sorted.bam
    samtools index -@ $t ${s}_nanopore_sorted.bam
else
    echo "Error: BAM files ${s}_illumina_sorted.bam and/or ${s}_nanopore_sorted.bam were not found." >&2
    exit
fi

# Merge BAM file ===============
# Sort SAM files by read names rather than mapped positions to solve the error of "No @HD tag found".
echo 'Merging BAM files'
samtools view -S -b -u -@ $t ${s}_illumina.sam | samtools sort -n --no-PG -@ $t --output-fmt BAM -o ${s}_illumina_name_sorted.bam -
samtools view -S -b -u -@ $t ${s}_nanopore.sam | samtools sort -n --no-PG -@ $t --output-fmt BAM -o ${s}_nanopore_name_sorted.bam -
if [ -f "${s}_illumina_name_sorted.bam" ] && [ -f "${s}_nanopore_name_sorted.bam" ]; then
    samtools merge -n -f --no-PG --output-fmt BAM -@ $t --reference $r -u -o ${s}_hybrid.bam ${s}_illumina_name_sorted.bam ${s}_nanopore_name_sorted.bam  # It's unlikely that BAM files are not produced in this and next commands at this stage.
    samtools sort -@ $t --output-fmt BAM -o ${s}_hybrid_sorted.bam ${s}_hybrid.bam
    samtools index -@ $t ${s}_hybrid_sorted.bam
    rm ${s}_illumina.sam ${s}_nanopore.sam ${s}_illumina_name_sorted.bam ${s}_nanopore_name_sorted.bam ${s}_hybrid.bam  # Clean-up
    echo 'Successfully went through this pipeline.'
else
    echo "Error: BAM files ${s}_illumina_name_sorted.bam and/or ${s}_nanopore_name_sorted.bam were not found." >&2
fi
