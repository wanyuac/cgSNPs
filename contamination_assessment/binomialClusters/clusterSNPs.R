#!/usr/bin/env Rscript
# Inference of multiplicity of bacterial isolates in a sample using unsupervised machine learning.
# Commandline: Rscript clusterSNPs.R --het [hetSNP table] --hom [homSNP table] --delim "\t" --eps [cut-off for posterior probabilities] --sampleName [sample name] --outdir [directory path not ended with a forward slash] --suffix SNPclusters
# Example command:
#   Rscript clusterSNPs.R --het hetSNP__info.tsv --hom homSNP__info.tsv --sampleName ST258 --eps 0.01 --outdir qc
# 
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# Edition: 3 May 2017

##### Define functions ###############
library(optparse)
library(flexmix)

optionList <- list(
    make_option("--hetSNP", type = "character", default = NULL, help = "A delimited text file of heterozygous SNPs", metavar = "character"),
    make_option("--homSNP", type = "character", default = NULL, help = "A delimited text file of homozygous SNPs [default= %default]", metavar = "character"),
    make_option("--delim", type = "character", default = "\t", help = "Deliminator of the input SNP tables [default= %default]", metavar = "character"),
    make_option("--eps", type = "character", default = "1e-4", help = "Minimum posterior probability below which probabilities are considered as zero [default= %default]", metavar = "character"),
    make_option("--sampleName", type = "character", default = "sample", help = "Sample name [default = %default]", metavar = "character"),
    make_option("--outdir", type = "character", default = ".", help = "Output directory (not ended with a forward slash) for plots [default = %default]", metavar = "character"),
    make_option("--suffix", type = "character", default = "SNPclusters", help = "Suffix of the output file name [default = %default]", metavar = "character")
)

optParser <- OptionParser(option_list = optionList);

parseDP4 <- function(x) {
    # parse the DP4 field in the SNP table. This is a subordinate function of mkCountMatrix.
    n <- nrow(x)
    x <- cbind(x, data.frame(DP.REF = rep(0, times = n), DP.ALT = rep(0, times = n)))
    
    for (i in 1 : n) {
        depths <- as.integer(unlist(strsplit(x$DP4[i], split = ",", fixed = TRUE)))
        x[i, c("DP.REF", "DP.ALT")] <- c(depths[1] + depths[2], depths[3] + depths[4])
    }
    
    return(x)
}

mkCountMatrix <- function(het, hom, delim = "\t") {
    # construct a two column matrix for classification of SNPs based on admixed binomial distributions using the flexmix package
    # first column: number of successes (reads showing an alternative allele); second column: number of failures (reads showing the reference allele)
    y <- parseDP4(read.delim(het, sep = delim, stringsAsFactors = FALSE)[, c("POS", "REF", "ALT", "DP", "DP4")])[, c("DP.ALT", "DP.REF")]
    #y <- rbind(y, parseDP4(read.delim(hom, sep = delim, stringsAsFactors = FALSE)[, c("POS", "REF", "ALT", "DP", "DP4")])[, c("DP.ALT", "DP.REF")])
    
    return(as.matrix(y))
}

# The following function is adapted from the script binomial_mixture.R in the MOIMIX package (https://github.com/bahlolab/moimix).
# Fit binomial mixture model on coverage data.
# Use functions initFlexmix and FLXMRglm in the flexmix package.
# Arguments:
#   y: a two-column matrix of read depths of the alt and ref alleles (no indels), respectively, without a header
#   k: a vector of mixture components to fit
#   n.iter: maximal number of iterations to run before getting parameter estimates converged
#   n.rounds: number of independent runs launched for each k to overcome local optimisation for the E-M algorithm.
#   min.prior: minimal prior probability below which the relevant component will be discarded for model fitting.
binommix <- function(y, k = 1 : 4, n.rounds = 5, n.iter = 1000, min.prior = 0) {
    # The first two arguments "y ~ 1" and "model = ..." will be passed directly to the function flexmix.
    # So do not change their positions in the argument list of stepFlexmix.
    fits <- stepFlexmix(y ~ 1, model = FLXMRglm(y ~ 1, family = "binomial"),
                        k = k, control = list(iter.max = n.iter, minprior = min.prior),
                        nrep = n.rounds)  # minprior = 0 turns off component removal
    
    return(fits)
}

summariseBinommix <- function(md, eps = 1e-4, sample.name) {  # md: the model fitted using binommix
    # The total number of components under all k's determines the final number of rows of m.
    m <- NULL
    if (class(md) == "stepFlexmix") {  # multiple k's
        ks <- names(md@models)  # number of components (k) tested in the function binommix; variable class: character
        for (k in ks) {
            s <- summary(md@models[[k]], eps = eps)
            names(s@comptab)[3] <- "post"  # comptab is a data frame. Replace "post>0" with "post".
            n <- nrow(s@comptab)  # number of components under the current k
            rownames(s@comptab) <- NULL  # They were "Comp.1", "Comp.2", etc.
            m.new <- cbind(sample = rep(sample.name, times = n), k = rep(k, times = n),
                           comp = 1 : n, s@comptab, BIC = rep(s@BIC, times = n),
                           AIC = rep(s@AIC, times = n), logLik = rep(s@logLik, times = n),
                           df = rep(md@models[[k]]@df, times = n), stringsAsFactors = FALSE)  # embed comptab into a larger data frame
            m <- rbind(m, m.new)
        }
    } else {  # class(md) = flexmix when a user only runs binommix for a single k
        s <- summary(md)
        names(s@comptab)[3] <- "post"
        n <- nrow(s@comptab)
        rownames(s@comptab) <- NULL
        m <- data.frame(sample = rep(sample.name, times = n), k = rep(md@k, times = n),
                        comp = 1 : n, s@comptab, BIC = rep(s@BIC, times = n),
                        AIC = rep(s@AIC, times = n), logLik = rep(s@logLik, times = n),
                        df = rep(md@df, times = n), stringsAsFactors = FALSE)
    }
    
    return(m)
}

##### Main program ###############
opts <- parse_args(optParser)
if (is.null(opts$hetSNP) | is.null(opts$homSNP)) {
    stop("Input files for heterozygous and/or homozygous SNPs are missing.")
}
eps <- as.numeric(opts$eps)  # By default, flexmix::summary uses 1e-4 as the cut-off for posterior probabilities.
outdir <- opts$outdir
if ((outdir != ".") & (!file.exists(outdir))) {
    dir.create(outdir)
}

# import read depths of the alternative allele(s) and the reference allele for a sample
md <- binommix(y = mkCountMatrix(het = opts$hetSNP, hom = opts$homSNP, delim = opts$delim))  # Under the default function setting, it fits models under k = 1..4 and launch 5 iterations for each k (20 runs in total).
md.stats <- summariseBinommix(md = md, eps = eps, sample.name = opts$sampleName)

# print results
print(paste0("Header: ", paste(names(md.stats), collapse = ",")))
write.table(md.stats, file = paste(outdir, paste0(opts$sampleName, "__", opts$suffix, ".tsv"), sep = "/"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)  # not printing the header as it will be easier to concatenate outputs of different samples
