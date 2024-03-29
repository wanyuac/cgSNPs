#!/usr/bin/env python

"""
This script is a Python3-based, partial reimplementation of my R script mergeGenomicRegions.R
in the repository BINF_toolkit (github.com/wanyuac/BINF_toolkit). It aims to exclude identical
regions defined in the output of script filterCoords.py and merge overlapping regions by their
coordinates. Compared to mergeGenomicRegions.R, this script does not look for complementary
regions after region merge, which makes its functionality more focused. User may use the output
coordinate list and bcftools to filter out variants identified in certain regions.

Input: A two-column, comma delimited (CSV), no-header text file in which each row defines the
start and end coordinates of a genomic region, respectively. The input can also come from the stdin.

Output format depends on parameter -f:
    csv, a test file of the same format as the input, but has its duplicated rows excluded and overlapping regions merged.
         Parameter '-n' is not used when -f = 'csv'.
    tsv, Conventional tab-delimited, 1-based coordinates (Chr\tStart\tEnd), for bcftools to exclude certain regions from the
         read alignment. Cf., option '--regions-file' of bcftools (http://samtools.github.io/bcftools/bcftools.html#common_options).
    bed, print coordinates in BED format (0-based) rather than 1-based position format in the output file.

Example command (three ways to run this script):
    python mergeGenomicRegions.py -i repeats.coords -o repeats_merged.coords
    python filterCoords.py -i genome.coords -I 90 | python mergeGenomicRegion.py -f csv -o repeats_merged.coords
    python filterCoords.py -i genome.coords -I 90 | python mergeGenomicRegion.py -f bed -n NZ_HG326223.1 -o repeats_merged.bed
    python filterCoords.py -i genome.coords -o repeats.coords -I 90 && python mergeGenomicRegion.py -i repeats.coords -f tsv -n NZ_HG326223.1 -o repeats_merged.tsv

Dependencies: Python 3 and packages pandas (pandas.pydata.org/)
Copyright 2020-2024 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 25 Apr 2020; last modification: 5 Feb 2024.
"""

import sys
import pandas as pd
import itertools as itr
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description = "Exclude duplicated region definitions and merge overlapping regions")
    parser.add_argument("-i", type = str, required = False, default = "", help = "Path to the input CSV file")
    parser.add_argument("-o", type = str, required = True, help = "Output file")
    parser.add_argument("-f", type = str, required = False, default = "tsv", help = "Output format: tsv (default)/bed/csv")
    parser.add_argument("-n", type = str, required = False, default = "Seq", help = "Sequence name (often the accession number) for the BED-format output")

    return parser.parse_args()


def main():
    args = parse_arguments()
    
    # Import the raw coordinate table
    try:
        if args.i == "":
            tab = pd.read_csv(sys.stdin, header = None, names = ["From", "To"], index_col = None, dtype = int)
        else:
            tab = pd.read_csv(args.i, header = None, names = ["From", "To"], index_col = None, dtype = int)
    except:
        print("The input file is not accessible.")

    # Remove duplicated rows: equivalent to unique(tab) in R.
    tab.drop_duplicates(subset = ["From", "To"], keep = "first", inplace = True)  # len(tab.index) shows a reduction in the row number.

    # Sort regions in an ascending order
    tab.sort_values(by = ["From", "To"], ascending = True, inplace = True)
    tab.reset_index(drop = True, inplace = True)

    # Merge overlapping regions
    n = len(tab.index)  # Number of remaining rows (regions) after deduplication
    if n > 1:
        tab_merge = mergeRegions(tab, n)
    else:
        tab_merge = tab  # Do not need to merge regions when there is one region left.

    # Prepare the output under a given format
    if args.f == "csv":
        # The simplest scenario: write the merged table in a CSV file
        tab_merge.to_csv(args.o, header = False, index = False, sep = ",")
    else:
        tab_merge.reset_index(drop = True, inplace = True)  # Discard indices for the join function
        chr = pd.DataFrame({"Seq" : list(itr.repeat(args.n, tab_merge["From"].count()))}, dtype = str)  # Create a new column
        tab_merge = chr.join(tab_merge)  # Equivalently append the new column to tab_merge    
        
        if args.f == "bed":
            """
            Some software requires the BED format (http://genome.ucsc.edu/FAQ/FAQformat#format1) for input coordinates. The BED format
            adopts a '0-start, half-open' format rather than the intuitive 'position' format for the ease of calculating region lengths.
            A detailed comparison of both formats can be found in a blog post:
            http://genome.ucsc.edu/blog/the-ucsc-genome-browser-coordinate-counting-systems/.
            The following command converts the intuitive '1-start, fully-closed' position format to the BED format. It substracts one
            from start positions and keeps end positions the same, as introduced in the blog post.
            """
            tab_merge["From"] -= 1  # This is an essential difference between BED files and TSV files.
        
        tab_merge.to_csv(args.o, header = False, index = False, sep = "\t", float_format = None)
    return


def mergeRegions(tab, nr):
    """
    This function merges regions [x, y1] and [y2, z] into [x, z] if y2 <= y1, assuming y1 > x, z > y2,
    and z > x.
    """
    n = 0  # Number of rows appended to the output data frame
    r1 = tab.iloc[n]  # Extract the first row (region) from the data frame and get a Pandas Series object
    ub1 = r1["To"]  # Upper boundary of the current region
    i_max = nr - 1  # Index of the last region
    for i in range(1, nr):
        r2 = tab.iloc[i]  # Get boundaries of the second region
        lb2 = r2["From"]  # Lower boundary of the new region
        ub2 = r2["To"]  # Upper boundary of the second region
        
        # Check validity
        if lb2 > ub2:
            sys.exit("Error: lower bound %d > upper bound %d." % (lb2, ub2))
        
        # Decision: either merge two regions or add a separate region.
        if lb2 <= (ub1 + 1):
            """
            Two regions overlap or immediately adjacent: merge them into a single one. Here, lb1 <= lb2
            because the data frame "tab" has been sorted in an ascending order.
            """
            if ub2 > ub1:
                """
                Extend the previous region r1 to ub2; else, do nothing as r2 is embedded in r1.
                """
                ub1 = ub2
                r1["To"] = ub1  # Extend the current region r1 and search for the next region
            
            if i == i_max:  # The last region is reached.
                if n == 0:  # Has not written any row yet
                    tab_merg = pd.DataFrame(data = {"From": [r1["From"]], "To": [r1["To"]]})  # Initialise the output data frame
                else:
                    tab_merg = pd.concat([tab_merg, pd.DataFrame(r1).transpose()], ignore_index = True)
        else:  # No overlap
            if n == 0:
                """
                This is the first region to be saved in the output data frame: push region "r1" (may
                have been extended) into the stack of regions.
                """
                tab_merg = pd.DataFrame(data = {"From": [r1["From"]], "To": [r1["To"]]})  # Initialise the output data frame
            else:  # No overlap and n > 0
                tab_merg = pd.concat([tab_merg, pd.DataFrame(r1).transpose()], ignore_index = True)
            n += 1
            
            # Deal with the last region
            if i == i_max:  # n must > 0
                tab_merg = pd.concat([tab_merg, pd.DataFrame(r2).transpose()], ignore_index = True)
            else:
                r1 = r2
                ub1 = ub2  # Move the upper boundary to that of r2 for subsequent comparisons

    return tab_merg


if __name__ == "__main__":
    main()
