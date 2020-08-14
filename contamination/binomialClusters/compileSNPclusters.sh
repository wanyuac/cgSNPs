#!/bin/bash
# Concatenate tables of SNP cluters that are produced by clusterSNPs_slurm.sh into a single file.
# Usage: bash compileSNPclusters.sh [individual input files] [output file]
# Example: bash compileSNPclusters.sh snp_clusters/*.tsv > snp_clusters_all.tsv
#
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# Editions: 4 May 2017

printf "sample\tk\tcomp\tprior\tsize\tpost\tratio\tBIC\tAIC\tlogLik\tdf\n"
for f in "$@"; do
    cat $f
done
