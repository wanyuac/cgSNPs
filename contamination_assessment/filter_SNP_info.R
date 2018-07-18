# Filter out SNPs in a VCF file given a list of genomic regions.
# Example:
# Rscript filterSNPs.R 5150_1#2__info.tsv \
# /vlsci/SG0006/shared/wan/Kp/snps/ref/NTUH-K2044__AP006725__rm.coords \
# /vlsci/SG0006/shared/wan/Kp/snps/vcf_parsed/het/filtered/5150_1#2__info 5150_1#2
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# First edition: 22/4/2017; the latest edition: 22/3/2018

args <- commandArgs(trailingOnly = TRUE)
f.snps <- args[1]
f.regions <- args[2]
out <- args[3]  # path and prefix of outputs
sample <- args[4]  # must not be NULL

snps <- read.delim(f.snps, stringsAsFactors = FALSE)
if (nrow(snps) > 0) {
    is.indel <- snps$INDEL == "1"
    snps.excl <- snps[is.indel, ]
    snps <- snps[!is.indel, ]  # remove indel information (may be none)
    snps$pass <- rep(TRUE, times = nrow(snps))

    # Import the filter for SNP coordinates (must not be an empty file)
    excl.regions <- read.csv(f.regions, header = FALSE)
    names(excl.regions) <- c("Start", "End")

    # Filter SNPs by their coordinates
    for (i in 1 : nrow(snps)) {
        pos <- snps$POS[i]
        for (j in 1 : nrow(excl.regions)) {
            line <- excl.regions[j, ]
            if ((pos >= line$Start) & (pos <= line$End)) {
                snps$pass[i] <- FALSE
                break
            }
        }
    }
    snps.pass <- subset(snps, pass)[, -ncol(snps)]  # drop the last column "pass"
    snps.excl <- rbind.data.frame(snps.excl, subset(snps, !pass)[, -ncol(snps)])  # include indels; valid even when both data frames have no rows.
} else {  # Sometimes there is no SNP (heterozygous or homozygous) in a sample.
    print(paste0("Warning: no SNP is present in the input file from the sample ", sample, "."))
    snps.pass <- snps
    snps.excl <- snps  # Both snps.pass and snps.excl are of zero row.
}

# Save results
if (nrow(snps.pass) > 0) {
    write.table(snps.pass, file = paste0(out, "_retained.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    n.pass <- nrow(snps.pass)
} else {
    print(paste("No heterozygous/homozygous SNPs in the sample", sample, "passed the filter.", sep = " "))
    n.pass <- 0
}

if (nrow(snps.excl) > 0) {
    write.table(snps.excl, file = paste0(out, "_removed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
    print(paste("No heterozygous/homozygous SNPs in the sample", sample, "got excluded by the filter.", sep = " "))
}

write.table(data.frame(sample = sample, snp.n = n.pass), file = paste0(out, "_snpNum.tsv"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)  # The number may be zero.
