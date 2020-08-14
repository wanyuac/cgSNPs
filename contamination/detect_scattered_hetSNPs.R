# Detecting scattered (unclustered) hetozygous SNPs (hetSNPs) using binomial tests
# Copyright 2018 Yu Wan <wanyuac@gmail.com>
# Licensed under the Apache License, Version 2.0
# Release: 19/6/2018, lastest edition: 25/6/2018

# DEFINE PARAMETERS ###############
library(optparse)

optionList <- list(
    make_option("--samples", type = "character", default = NULL, action = "store",
               help = "A text file containing a list of sample names", metavar = "character"),
    make_option("--indir", type = "character", default = "Het", action = "store",
                help = "Directory of parsed VCF files of SNPs. No end forward slash.",
                metavar = "character"),
    make_option("--suffix", type = "character", default = "__info_retained.tsv",
                action = "store", help = "Suffix of VCF files. For example, \'__info.tsv\' for \'NJST258__info.tsv\'",
                metavar = "character"),
    make_option("--genomeLen", type = "integer", default = 0, action = "store",
                help = "Length of the reference genome [default = %default]", metavar = "character"),
    make_option("--rmRegions", type = "character", default = NULL, action = "store",
                help = "(Optional) A two-column CSV file defining regions to be excluded from the SNP table"),
    make_option("--windowSizes", type = "character", default = "100,500,1000,1500",
                action = "store", help = "(Optional) Number of bases (comma delimited) in each sliding window [default = %default]",
                metavar = "character"),
    make_option("--maxP", type = "double", default = 0.05, action = "store",
                help = "Maximal Bonferroni-corrected p-value for significant SNP clusters [default = %default]",
                metavar = "character"),
    make_option("--cores", type = "integer", default = 1, action = "store",
                help = "Number of computational cores to run tasks [default = %default]",
                metavar = "character"),
    make_option("--outdir", type = "character", default = ".", action = "store",
                help = "(Optional) Output directory (not ended with a forward slash) [default = %default]",
                metavar = "character"),
    make_option("--outPrefix", type = "character", default = NULL, action = "store",
                help = "(Optional) Prefix for the output file [default = %default]",
                metavar = "character"),
    make_option("--override", type = "logical", action = "store_true", default = FALSE,
                help = "(Optional) Override saved previous results to launch binomial tests.",
                metavar = "logical"),
    make_option("--rmTemp", type = "logical", action = "store_true", default = FALSE,
                help = "(Optional) Remove temporary files at the final stage",
                metavar = "logical")
)

# DEFINE FUNCTIONS ###############
# Environmental setting ================
checkDir <- function(d) {
    if (!file.exists(d)) {  # check the presence of \temp
        dir.create(d)
    }
}

# Data preparation ===============
mergeRegions <- function(rs) {
    # Merge overlapping regions defined in a data frame.
    # First, get unique rows to reduce the number of comparisons.
    # This function is a modification of my script mergeGenomicRegions.R in the BINF_toolkit repository on GitHub.
    rs <- unique(rs)

    if (nrow(rs) > 1) {
        # Sort beginnings of regions so that their heads do not go backwards
        rs <- rs[order(rs$From, decreasing = FALSE), ]

        # Then, merge overlapping regions
        i <- 1  # initial row pointer for the output data frame
        r1 <- rs[i, ]  # initialise the output using the region closest to the 5' end as a beginning
        ub1 <- r1$To[i]  # upper bound of the current region

        for (j in 2 : nrow(rs)) {  # must guarantee there are >= 2 rows in the table
            r2 <- rs[j, ]  # take another row from the downstream side
            lb2 <- r2$From  # lower bound of the second region
            ub2 <- r2$To

            # There are only two behavious: either merge two regions or adding a separate region.
            if (lb2 <= (ub1 + 1)) { # When two regions are overlapping or adjacent: concatenate them into a single one.
                # Herein lb1 <= lb2 because the data frame has been sorted in an ascending order of "From".
                if (ub2 > ub1) {
                    # extend region 1 into region 2
                    r1$To[i] <- ub2
                    ub1 <- ub2  # move the upper bond towards the 3' end while the lb1 is kept the same
                }  # else, do nothing as the second range is a subset of the first one
            } else { # When the second resion is non-overlapping: push it into the stack of regions and refresh r1
                r1 <- rbind.data.frame(r1, r2)
                i <- i + 1  # move the focus to the next row in r1
                ub1 <- ub2  # lb1 is not useful.
            }
        }
    } else {  # Do nothing but copy rs to r1 when the deduplicated rs contains only a single row.
        r1 <- rs
    }

    r1$Bases <- r1$To - r1$From + 1  # number of bases in each region

    return(r1)
}

# Binomial tests ===============
parseDP4 <- function(x) {
    # Parse the DP4 field in the SNP table. Two new columns will be appended to the original data frame x.
    # This function returns fold coverage of alternative and reference bases by high-quality base calls.
    n <- nrow(x)
    x <- cbind.data.frame(x, data.frame(DP_REF = rep(0, times = n), DP_ALT = rep(0, times = n)))

    for (i in 1 : n) {
        depths <- as.integer(unlist(strsplit(x$DP4[i], split = ",", fixed = TRUE)))
        x[i, c("DP_REF", "DP_ALT")] <- c(depths[1] + depths[2], depths[3] + depths[4])
    }

    return(x)
}

testSNPdensities <- function(s, L0, ws, info_files, ci, tmp_dir, prefix, override) {
    # Assuming SNPs have been filtered with the exclusive regions already.
    info <- read.delim(info_files[[s]], stringsAsFactors = FALSE)  # import the SNP information
    info <- info[, c("POS", "DP4")]  # only use high-quality bases
    info <- info[order(info$POS, decreasing = FALSE), ]  # sort rows by positions in an ascending order
    info <- parseDP4(info)  # append two columns (DP_REF and DP_ALT) to this data frame
    info <- subset(info, (DP_REF > 0) & (DP_ALT > 0))  # remove hetSNPs caused by low-quality bases, whose B-allele frequency = 0 or 100% according to the DP4 information

    # Initialise the output
    tests <- vector(mode = "list", length = (4 + length(ws)))  # initialise the result list
    names(tests) <- c("sample", "snp_num", "info", "p", as.character(ws))
    tests[["sample"]] <- s  # The parallel function does not return a named list.
    tests[["info"]] <- info
    snp_num <- nrow(info)  # high-quality base calls
    tests[["snp_num"]] <- snp_num  # number of remaining SNPs

    if (snp_num > 0) {
        # Estimate probability of hetSNPs ("success") in the sum of remaining genomic regions
        p0 <- nrow(info) / L0

        # Perform binomial tests under each window size
        pos_min <- min(info$POS)
        pos_max <- max(info$POS)

        for (w in ws) {
            # Initialise the output file
            tmp <- paste(tmp_dir, paste(prefix, s, w, "binom.tsv", sep = "__"),
                         sep = "/")  # tmp_dir/prefix__s__w__binom.tsv
            if ((!file.exists(tmp) || override)) {
                w_start <- pos_min  # start position of the first window
                w_start_max <- pos_max - w + 1  # the begging position of the last window
                #t <- data.frame(ID = integer(0), From = integer(0), To = integer(0),
                #                Count = integer(0), P = numeric(0))

                # Since rbind.data.frame is slow when there are lots of rows
                # to combine, it may become the bottleneck of speed. Hence I
                # make the function to print the lines and read them as a
                # data frame afterwards.

                f <- file(tmp, "wt")  # make a connection to the output file to increase the accessing speed
                write(paste("ID", "From", "To", "Count", "P", sep = "\t"), file = f)  # the header line

                # Slide the window
                i <- 1  # window ID
                while(w_start <= w_start_max) {
                    w_end <- w_start + w - 1
                    n <- sum((info$POS >= w_start) & (info$POS <= w_end))  # number of hetSNPs within the current window
                    if (n > 0) {  # There is no need to perform the test as p < p0 always hold in this window.
                        # Null hypothesis: the true probability of observing a SNP at
                        # each base within this window is no more than p0.
                        # Since p > p0 when the window shows an elevated density of
                        # SNPs, we reject the null hypothesis when the p-value <= maxP.
                        # Also, the binom.test function always return a p-value of one
                        # when x = 0 regardless of hypotheses or p0 because the estimated
                        # probability of success equals zero.
                        p <- binom.test(x = n, n = w, p = p0, alternative = "greater",
                                        conf.level = ci)
                        #t <- rbind.data.frame(t, data.frame(ID = i, From = w_start,
                        #                                    To = w_end, Count = n,
                        #                                    P = p$p.value))
                        write(paste(i, w_start, w_end, n, round(p$p.value, digits = 8),
                                    sep = "\t"),
                              file = f, append = TRUE)
                    }
                    i <- i + 1
                    w_start <- w_start + 1  # move to the next window
                }
                close(f)
            } else {
                print(paste0("Reading previous results from ", tmp, "."))  # requires makeCluster(outfile = ...)
            }

            # Read the temporary file
            t <- read.delim(file = tmp)

            # Bonferroni correction for p-values
            t$P_adj <- p.adjust(p = t$P, method = "bonferroni")
            tests[[as.character(w)]] <- t
        }
        tests[["p"]] <- p0
    } else {  # when no hetSNPs left
        tests[["p"]] <- 0
    }

    # Return a list per sample
    return(tests)
}

# Summarise results ===============
reorganiseTestsList <- function(tests) {
    # Reorganise elements of the list tests and return a named list.
    n <- length(tests)
    out <- list()
    for (lst in tests) {
        out[[lst$sample]] <- lst[which(names(lst) != "sample")]
    }

    return(out)
}

countSNPsInRegions <- function(regions, snps) {
    # This is a subordinate function of getHighDensityRegions.
    # snps: an integer vector for SNP positions
    n <- nrow(regions)
    regions$SNP_num <- integer(length = n)
    for (i in 1 : n) {
        r <- regions[i, ]
        snps_in <- snps[which((snps >= r$From) & (snps <= r$To))]
        # Shrink the high-density region
        # It is impossible to have an empty snps_in vector in each region.
        regions[i, c("From", "To", "SNP_num")] <- c(min(snps_in), max(snps_in),
                                                    length(snps_in))
    }

    return(regions)
}

getHighDensityRegions <- function(tests, p_max) {
    # Concatenate adjacent windows of significantly higher SNP density into larger regions and count SNPs within each region
    n <- length(tests)
    out <- vector(mode = "list", length = n)
    samples <- names(tests)
    names(out) <- samples  # copy sample names as element names of the output list
    for (s in samples) {  # go through each sample
        t <- tests[[s]]  # The element t is a named list.
        snps <- t$info  # already filtered out hetSNPs caused by low base quality
        snps <- sort(snps$POS, decreasing = FALSE)  # an integer vector
        ws <- names(t)
        ws <- ws[-which(ws %in% c("snp_num", "info", "p"))]  # The element "sample" has been removed by the previous function reorganiseTestsList.
        if (t$snp_num > 0) {  # When t$snp_num = nrow(t$info) > 0, other elements cannot be NULL or empty data frames.
            out_s <- data.frame(W_size = integer(0), From = integer(0), To = integer(0),
                                Bases = integer(0), SNP_num = integer(0))  # the output for the current sample
            for (w in ws) {
                tw <- t[[w]]
                tw_sig <- subset(tw, P_adj <= p_max)
                if (nrow(tw_sig) > 0) {  # When there are high-density regions
                    tw_sig <- mergeRegions(tw_sig[, c("From", "To")])
                    tw_sig <- cbind.data.frame(W_size = rep(as.integer(w), times = nrow(tw_sig)),
                                               tw_sig)  # append a column of window size
                    tw_sig <- countSNPsInRegions(tw_sig, snps)
                    out_s <- rbind.data.frame(out_s, tw_sig)
                } else {  # when there is no high-density region
                    out_s <- rbind.data.frame(out_s,
                                              data.frame(W_size = as.integer(w),
                                                         From = 0, To = 0,
                                                         Bases = 0, SNP_num = 0))
                }
            }
        } else {
            # When all hetSNPs are caused by low-quality base calls, consider
            # the whole genome as of "high density".
            nw <- length(ws)  # number of window sizes
            out_s <- data.frame(W_size = as.integer(ws), From = rep(min(snps), times = nw),
                                To = rep(max(snps), times = nw))
            out_s$Bases <- out_s$To - out_s$From + 1
            out_s$SNP_num <- rep(length(snps), times = nw)  # All hetSNPs are treated as normal. Hence there is no scattered hetSNPs left.
        }
        out[[s]] <- out_s
    }

    return(out)
}

summariseSNPsPerSample <- function(spikes, tests) {
    # spikes: SNP clusters
    samples <- names(spikes)
    out <- data.frame(Sample = character(0), SNP_total = integer(0),
                      W_size = integer(0), P0 = numeric(0),
                      Clstr_num = integer(0), SNP_clstr = integer(0),
                      SNP_scat = integer(0), stringsAsFactors = FALSE)
    for (s in samples) {
        t <- tests[[s]]
        spk <- spikes[[s]]
        snp_num <- t$snp_num
        ws <- unique(spk$W_size)
        for (w in ws) {
            tw <- subset(spk, W_size == w)
            snp_clustered <- sum(tw$SNP_num)
            out <- rbind.data.frame(out, data.frame(Sample = s, SNP_total = snp_num,
                                                    W_size = w, P0 = t$p,
                                                    Clstr_num = nrow(tw),
                                                    SNP_clstr = snp_clustered,
                                                    SNP_scat = snp_num - snp_clustered,
                                                    stringsAsFactors = FALSE))
        }
    }

    return(out)
}

# READ AND PARSE ARGUMENTS ###############
# Read arguments
optParser <- OptionParser(option_list = optionList)
opts <- parse_args(optParser)

# Check if the genome size is valid
if (opts$genomeLen <= 0) {
    stop("Argument error: length of the reference genome must be positive.")
}

# Create output directories
tmp_dir <- paste(opts$outdir, paste0(opts$outPrefix, "_temp"), sep = "/")  # e.g., outdir/prefix_temp
checkDir(opts$outdir)
checkDir(tmp_dir)

# Read sample names
samples <- read.delim(file = opts$samples, header = FALSE, stringsAsFactors = FALSE)  # There is only a single column without a header in the file.
samples <- samples[, 1]

# Parse window sizes
ws <- as.integer(unlist(strsplit(x = opts$windowSizes, split = ",", fixed = TRUE)))

# Organise input files
info_files <- sapply(samples, function(s) paste(opts$indir, paste0(s, opts$suffix), sep = "/"))  # returns a named vector of filenames

# Confidence level
ci <- 1 - opts$maxP

# Import and merge overlapping exclusive regions
excl <- read.csv(file = opts$rmRegions, header = FALSE, stringsAsFactors = FALSE)
names(excl) <- c("From", "To")
excl <- mergeRegions(excl)  # Returns a data frame of three columns: "From", "To" and "Bases".
excl_bp <- sum(excl$Bases)  # number of bases to be excluded
L0 <- opts$genomeLen - excl_bp  # remaining bases for which hetSNPs are to be analysed
print(paste(paste0("[", Sys.time(), "]:"), "There are", excl_bp,
            "bp in the reference genome are masked from this assessment, leaving",
            L0, "bp for analsis.", sep = " "))

# PROCESS EACH SAMPLE UNDER EVERY WINDOW SIZE ###############
library(parallel)

# Determine number of cores
n_cores <- detectCores()
if (opts$cores < n_cores & opts$cores > 0) {
    n_cores <- opts$cores
}
print(paste(paste0("[", Sys.time(), "]:"),
            "Launching", n_cores, "cores to identify elevated SNP density for",
            length(samples), "samples.", sep = " "))

# Cluster identification for each sample in parallel
cl <- makeCluster(n_cores, outfile = paste(opts$outdir, paste(opts$outPrefix, "parallel.log",
                                                              sep = "__"), sep = "/"))
clusterExport(cl = cl,
              varlist = list("L0", "ws", "info_files", "ci", "parseDP4",
                             "tmp_dir", "opts"),
              envir = environment())  # make variables accessible to different cores
tests <- parLapply(cl, samples, testSNPdensities, L0, ws, info_files, ci, tmp_dir,
                   opts$outPrefix, opts$override)
stopCluster(cl)
print(paste(paste0("[", Sys.time(), "]:"), "One-sided binomial tests are finished.",
            sep = " "))
tests <- reorganiseTestsList(tests)

# PRODUCE SUMMARY STATISTICS ###############
# Merge adjacent windows into larger regions
spikes <- getHighDensityRegions(tests = tests, p_max = opts$maxP)

# Summary statistics per sample
hetSNP_summary <- summariseSNPsPerSample(spikes, tests)

# SAVE RESULTS ###############
# Since the list "tests" is often extremely large, results are not saved in a
# single file.
print(paste0("[", Sys.time(), "]: Saving results."))

if (is.null(opts$outPrefix)) {
    saveRDS(excl, file = paste(opts$outdir, "excl.rds", sep = "/"))
    saveRDS(spikes, file = paste(opts$outdir, "clusters.rds", sep = "/"))
    saveRDS(hetSNP_summary, file = paste(opts$outdir, "summary.rds", sep = "/"))
    saveRDS(tests, file = paste(opts$outdir, "tests.rds", sep = "/"))  # The slowest step goes last so users can download other files in the mean time.
} else {
    saveRDS(excl, file = paste(opts$outdir,
                               paste(opts$outPrefix, "excl.rds", sep = "__"), sep = "/"))
    saveRDS(spikes, file = paste(opts$outdir,
                                 paste(opts$outPrefix, "clusters.rds", sep = "__"), sep = "/"))
    saveRDS(hetSNP_summary, file = paste(opts$outdir,
                                         paste(opts$outPrefix, "summary.rds", sep = "__"), sep = "/"))
    saveRDS(tests, file = paste(opts$outdir,
                                paste(opts$outPrefix, "tests.rds", sep = "__"), sep = "/"))
}

# CLEARANCE ###############
if (opts$rmTemp) {
    system(paste0("rm -rf ", tmp_dir), wait = TRUE)
}

print(paste0("[", Sys.time(), "]: Done."))
