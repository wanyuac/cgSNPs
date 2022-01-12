#!/usr/bin/env Rscript
# Generate a midpoint-rooted neighbour-joining (NJ) tree from a distance matrix.
#
# Input: a matrix of SNP distances from software snp-dists
# Output: a Newick-format NJ tree
# Example: Rscript --vanilla dist2tree.R --i snp.dist --o snp_dist.newick --r ref_genome
# Dependencies: R packages optparse, ape, phytools
#
# Copyright (C) 2022 Yu Wan <wanyuac@126.com>
# Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
# Creation: 11 Jan 2022; the latest update: 11 Jan 2022.

library(optparse)
library(ape)
library(phytools)

# Options ###############
options <- list(
    make_option("--i", type = "character", help = "A tab-delimited matrix from snp-dists or compatible software"),
    make_option("--o", type = "character", default = "snp_dist.newick", help = "Output tree file [default %default]"),
    make_option("--r", type = "character", default = "reference", help = "(Optional) Name for the reference genome [default %default]")
)

# Main ###############
opt.parser <- OptionParser(usage = "Rscript %prog [options]", option_list = options)
opts <- parse_args(opt.parser)

message(paste("Reading distance matrix.", opts$i, sep = " "))
D <- as.matrix(read.delim(file = opts$i, row.names = 1))
r <- opts$r
if (r != "reference") {
    j <- which(rownames(D) == "reference")
    k <- which(colnames(D) == "reference")
    if (length(j) == 0 || length(k) == 0 || j != k) {
        stop("Rows and colums in the distance matrix do not completely match. No tree will be generated.")
    } else {
        rownames(D)[j] <- r
        colnames(D)[k] <- r
    }
}
message(paste("Generating a neighbour-joining tree for", nrow(D), "taxa.", sep = " "))
write.tree(phy = ladderize(midpoint.root(nj(D)), right = FALSE), file = opts$o)
message(paste("Tree file", opts$o, "was generated.", sep = " "))