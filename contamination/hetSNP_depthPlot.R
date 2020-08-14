#!/usr/bin/env Rscript
# Draw read depths of heterozygous (het) and homozygous (hom) SNPs for a single sample.
# Read depths of heterozygous SNPs will be drawn as foreground points; those of homozygous SNPs will be drawn as background points.
# Commandline: Rscript hetSNP_depthPlot.R --het [hetSNP table] --hom [homSNP table] --delim "\t" --sampleName [sample name] --genomeLen [length (bp) of the reference genome] \
# --outdir [directory path not ended with a forward slash] --suffix [filename suffix]
# Example command:
#   Rscript drawHetSNPDP.R --het hetSNP__info.tsv --hom homSNP__info.tsv --sampleName ST258_1 --genomeLen 5.25e6 --outdir pic --suffix hetSNPs
# Notice because RedDog only store variants in VCFs, there is no band of BAF = 0 for homozygous SNPs in the BAF plot.
#
# Copyright 2017 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# Edition: 28/4, 2/5/2017, 8/3/2018; last modification: 13/8/2020

##### Arguments ###############
library(optparse)

optionList <- list(
    make_option("--hetSNP", type = "character", default = NULL, help = "A delimited text file of heterozygous SNPs", metavar = "character"),
    make_option("--homSNP", type = "character", default = NULL, help = "A delimited text file of homozygous SNPs [default= %default]", metavar = "character"),
    make_option("--delim", type = "character", default = "\t", help = "Deliminator of the input SNP tables [default= %default]", metavar = "character"),
    make_option("--sampleName", type = "character", default = "sample", help = "Sample name [default = %default]", metavar = "character"),
    make_option("--genomeLen", type = "integer", default = 5.25e6, help = "Length of the reference genome [default = %default]", metavar = "integer"),
    make_option("--outdir", type = "character", default = ".", help = "Output directory (not ended with a forward slash) for plots [default = %default]", metavar = "character"),
    make_option("--suffix", type = "character", default = "hetSNPs", help = "Suffix of the output file name [default = %default]", metavar = "character")
)

optParser <- OptionParser(option_list = optionList)


##### Constants ###############
IMG.WIDTH <- 1000
IMG.HEIGHT.SINGLE <- 600
IMG.HEIGHT.TRIPLE <- 1200
RES <- 72  # or 300 for journal figures
UNIT <- "px"  # or "mm
COL.REF <- "blue"  # colour for data points of reference alleles
COL.ALT <- "red"  # colour for data points of alternative alleles
COL.HOM <- "grey70"  # colour for data points of homozygous SNPs in the third panel
COL.HET <- "red"  # colour for data points represting heterozygous SNPs in the third panel
COL.AVGDP <- "blue"  # colour for the horizontal line of the average read depth in the third panel
POINT.SHAPE <- 16  # the "pch" argument for plots
POINT.SIZE <- 0.75  # the "cex" argument for points
X.INTERVAL <- 0.5e6  # width of every interval of the X axis (500 kb by default)
DP.INTERVAL <- 20
CEX.AXIS <- 1.5
CEX.AXIS.LABEL <- 2


##### Data structure and functions ###############
checkColumnNames <- function(x, target.names, id) {
    presence <- target.names %in% names(x)
    if (sum(presence) < length(target.names)) {
        stop(paste("Missing column name(s) in the SNP table", paste0(id, ":"), paste(target.names[!presence], collapse = ","), sep = " "))
    }
}

parseDP4 <- function(x) {
    # parse the DP4 field in the SNP table. Two new columns will be appended to the original data frame x.
    n <- nrow(x)
    x <- cbind(x, data.frame(DP.REF = rep(0, times = n), DP.ALT = rep(0, times = n)))

    for (i in 1 : n) {
        depths <- as.integer(unlist(strsplit(x$DP4[i], split = ",", fixed = TRUE)))
        x[i, c("DP.REF", "DP.ALT")] <- c(depths[1] + depths[2], depths[3] + depths[4])
    }

    return(x)
}


##### Import data ###############
opts <- parse_args(optParser)  # parse arguments

# check the presence of input tables
if (is.null(opts$hetSNP)) {
    stop("Argument error: at least a table of heterozygous SNPs must be provided.")
}

# check the output directory
if ((opts$outdir != ".") & (!file.exists(opts$outdir))) {
    dir.create(opts$outdir)
}

# import SNP tables
het <- read.delim(opts$hetSNP, sep = opts$delim, stringsAsFactors = FALSE)
names(het) <- toupper(names(het))
checkColumnNames(het, c("POS", "DP", "DP4"), "het")
het <- parseDP4(het)  # parse the DP4 column for read depths
n.het <- nrow(het)  # number of hetSNPs
het.maxDP <- ceiling(max(het$DP.REF, het$DP.ALT) / 10) * 10  # maximal read depth of hetSNPs

if (is.null(opts$homSNP)) {
    no.homSNP <- TRUE
} else {
    no.homSNP <- FALSE
    hom <- read.delim(opts$homSNP, sep = opts$delim, stringsAsFactors = FALSE)
    names(hom) <- toupper(names(hom))
    checkColumnNames(hom, c("POS", "DP", "DP4"), "hom")
    hom <- parseDP4(hom)
    het$DP <- het$DP.REF + het$DP.ALT  # I deliberately do not use the original het$DP for plots as it contains read depths of other alleles when the alternative allele is not unique for a SNP.
    hom$DP <- hom$DP.REF + hom$DP.ALT
    avgDP <- round(mean(c(hom$DP, het$DP)), digits = 2)  # average read depth of SNPs based on high-quality base calls
    snp.maxDP <- ceiling(max(c(het$DP, hom$DP)) / DP.INTERVAL) * DP.INTERVAL  # define the upper Y limit for the second panel

    # calculate B-allele (i.e. the alternative allele) frequencies
    het$BAF <- round(het$DP.ALT / het$DP, digits = 4)
    hom$BAF <- round(hom$DP.ALT / hom$DP, digits = 4)
}


###### Set up arguments of plots ###############
img.name <- paste(opts$outdir, paste(n.het, opts$sampleName, paste0(opts$suffix, ".png"), sep = "__"),
                  sep = "/")  # output: [outdir]/[number of hetSNPs]__[sample name]__[suffix].png

# configure ticks and labels of the X axis (SNP positions)
if (opts$genomeLen %% X.INTERVAL > 0) {  # for the most of cases
    x.ticks <- c(seq(0, opts$genomeLen %/% X.INTERVAL * X.INTERVAL, by = X.INTERVAL), opts$genomeLen)
    x.labels <- c("0", paste0(x.ticks[2 : (length(x.ticks) - 1)] / 1000, "k"), "")
} else {
    x.ticks <- seq(0, opts$genomeLen, by = X.INTERVAL)
    x.labels <- c("0", paste0(x.ticks[-1], "k"))
}


##### Plot one/two panels ###############
if (no.homSNP) {  # a single panel: hetSNPs
    png(filename = img.name, width = IMG.WIDTH, height = IMG.HEIGHT.SINGLE, res = RES, units = UNIT)
    par(mar = c(4.2, 4.2, 0.1, 0.1), las = 1)
    plot(x = het$POS, y = het$DP.REF, col = COL.REF, pch = POINT.SHAPE, cex = POINT.SIZE,
         xlim = c(0, opts$genomeLen), ylim = c(0, het.maxDP), axes = FALSE,
         xlab = "Genomic position (bp)", ylab = "Read depth per allele", cex.lab = CEX.AXIS.LABEL - 0.5)  # initialise the plotting area
    points(x = het$POS, y = het$DP.ALT, col = COL.ALT, pch = POINT.SHAPE, cex = POINT.SIZE)
    text(x = 1e5, y = het.maxDP, labels = paste0("n = ", n.het))  # print the number of hetSNPs at the top-left corner of this panel
    axis(side = 1, at = x.ticks, labels = x.labels, cex.axis = CEX.AXIS - 0.25)
    axis(side = 2, at = c(seq(0, het.maxDP, by = DP.INTERVAL)), cex.axis = CEX.AXIS - 0.25)
} else {  # three panels, for hetSNPs, BAF of all SNPs and read depth per SNP, respectively
    png(filename = img.name, width = IMG.WIDTH, height = IMG.HEIGHT.TRIPLE, res = RES, units = UNIT)
    layout(matrix(c(1, 2, 3), nrow = 3, ncol = 1), widths = c(1, 1, 1), heights = c(1, 1, 1))
    par(mar = c(4.2, 4.2, 0.1, 0.1), mgp = c(2.5, 0.8, 0), las = 1)

    # a. the top panel for absolute read depths of hetSNPs
    plot(x = het$POS, y = het$DP.REF, col = COL.REF, pch = POINT.SHAPE, cex = POINT.SIZE,
         xlim = c(0, opts$genomeLen), ylim = c(0, het.maxDP), axes = FALSE,
         xlab = "", ylab = "Read depth per allele", cex.lab = CEX.AXIS.LABEL)  # initialise the plotting area
    points(x = het$POS, y = het$DP.ALT, col = COL.ALT, pch = POINT.SHAPE, cex = POINT.SIZE)
    text(x = 1e5, y = het.maxDP, labels = paste0("n = ", n.het), cex = CEX.AXIS)
    axis(side = 1, at = x.ticks, labels = x.labels, cex.axis = CEX.AXIS)
    axis(side = 2, at = c(seq(0, het.maxDP, by = DP.INTERVAL)), cex.axis = CEX.AXIS, line = -1)

    # b. draw a BAF plot in the middle panel
    # B-allele frequency: Alternative-allele frequency
    plot(x = hom$POS, y = hom$BAF, col = COL.HOM, pch = POINT.SHAPE, cex = POINT.SIZE,
         xlim = c(0, opts$genomeLen), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "B-allele frequency", cex.lab = CEX.AXIS.LABEL)
    points(x = het$POS, y = het$BAF, col = COL.ALT, pch = POINT.SHAPE, cex = POINT.SIZE)
    axis(side = 1, at = x.ticks, labels = x.labels, cex.axis = CEX.AXIS)
    axis(side = 2, at = seq(0, 1, by = 0.2), cex.axis = CEX.AXIS, line = -1)

    # c. the bottom panel for comparing read depths of hetSNPs with those of homSNPs
    plot(x = hom$POS, y = hom$DP, col = COL.HOM, pch = POINT.SHAPE, cex = POINT.SIZE,
         xlim = c(0, opts$genomeLen), ylim = c(0, snp.maxDP), axes = FALSE,
         xlab = "Genomic position (bp)", ylab = "Read depth per SNP", cex.lab = CEX.AXIS.LABEL)
    abline(h = avgDP, col = COL.AVGDP, lty = 1)  # use a dashed line to denote the average read depth of SNPs
    points(x = het$POS, y = het$DP, col = COL.HET, pch = POINT.SHAPE, cex = POINT.SIZE)  # finally, overlay points of hetSNPs on the figure
    text(x = 3e5, y = snp.maxDP, labels = paste0("avg. depth = ", avgDP), cex = CEX.AXIS)  # print the average read depth of all SNPs at the top-left corner
    axis(side = 1, at = x.ticks, labels = x.labels, cex.axis = CEX.AXIS)
    axis(side = 2, at = c(seq(0, snp.maxDP, by = DP.INTERVAL)), cex.axis = CEX.AXIS, line = -1)
}

dev.off()
