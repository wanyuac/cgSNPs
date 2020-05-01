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

Output: A test file of the same format as the input, but has its duplicated rows excluded and
overlapping regions merged.

Example command (three ways to run this script):
    python mergeGenomicRegions.py -i repeats.coords -o repeats_merged.coords
    python filterCoords.py -i genome.coords -I 90 | python mergeGenomicRegion.py -o repeats_merged.coords
    python filterCoords.py -i genome.coords -o repeats.coords -I 90 && python mergeGenomicRegion.py -i repeats.coords -o repeats_merged.coords

Dependencies: Python v3 and packages pandas (pandas.pydata.org/)
Copyright 2020 Yu Wan <wanyuac@126.com>
Licensed under the GNU General Public Licence version 3 (GPLv3) <https://www.gnu.org/licenses/>.
Publication: 25 Apr 2020; last modification: 1 May 2020.
"""

import sys
import pandas as pd
from argparse import ArgumentParser

def parse_arguments():
    parser = ArgumentParser(description = 'Exclude duplicated region definitions and merge overlapping regions')
    parser.add_argument('-i', type = str, required = False, default = '', help = 'Path to the input CSV file')
    parser.add_argument('-o', type = str, required = True, help = 'Path to the output CSV file')

    return parser.parse_args()


def main():
    args = parse_arguments()
    
    # Import the raw coordinate table
    try:
        if args.i == '':
            tab = pd.read_csv(sys.stdin, header = None, names = ['From', 'To'], index_col = None, dtype = int)
        else:
            tab = pd.read_csv(args.i, header = None, names = ['From', 'To'], index_col = None, dtype = int)
    except:
        print('The input file is not accessible.')

    # Remove duplicated rows: equivalent to unique(tab) in R.
    tab.drop_duplicates(subset = ['From', 'To'], keep = 'first', inplace = True)  # len(tab.index) shows a reduction in the row number.

    # Sort regions in an ascending order
    tab.sort_values(by = ['From', 'To'], ascending = True, inplace = True)
    tab.reset_index(drop = True, inplace = True)

    # Merge overlapping regions
    n = len(tab.index)  # Number of remaining rows (regions)
    if n > 1:
        tab_merg = mergeRegions(tab, n)
    else:
        tab_merg = tab  # Do not need to merge regions when there is one region left.

    # Write the merged table in a CSV file
    tab_merg.to_csv(args.o, header = False, index = False)

    return


def mergeRegions(tab, nr):
    """
    This function merges regions [x, y1] and [y2, z] into [x, z] if y2 <= y1, assuming y1 > x, z > y2,
    and z > x.
    """
    n = 0  # Number of rows appended to the output data frame
    r1 = tab.iloc[n]  # Extract the first row (region) from the data frame and get a Pandas Series object
    ub1 = r1['To']  # Upper boundary of the current region
    i_max = nr - 1  # Index of the last region
    for i in range(1, nr):
        r2 = tab.iloc[i]  # Get boundaries of the second region
        lb2 = r2['From']  # Lower boundary of the new region
        ub2 = r2['To']  # Upper boundary of the second region
        
        # Check validity
        if lb2 > ub2:
            sys.exit('Error: lower bound %d > upper bound %d.' % (lb2, ub2))
        
        # Decision: either merge two regions or add a separate region.
        if lb2 <= (ub1 + 1):
            """
            Two regions overlap or immediately adjacent: merge them into a single one. Here, lb1 <= lb2
            because the data frame 'tab' has been sorted in an ascending order.
            """
            if ub2 > ub1:
                """
                Extend the previous region r1 to ub2; else, do nothing as r2 is embedded in r1.
                """
                ub1 = ub2
                r1['To'] = ub1  # Extend the current region r1 and search for the next region
            
            if i == i_max:  # The last region is reached.
                if n == 0:  # Has not written any row yet
                    tab_merg = pd.DataFrame(data = {'From': [r1['From']], 'To': [r1['To']]})  # Initialise the output data frame
                else:
                    tab_merg = tab_merg.append(r1)
        else:  # No overlap
            if n == 0:
                """
                This is the first region to be saved in the output data frame: push region 'r1' (may
                have been extended) into the stack of regions.
                """
                tab_merg = pd.DataFrame(data = {'From': [r1['From']], 'To': [r1['To']]})  # Initialise the output data frame
            else:  # No overlap and n > 0
                tab_merg = tab_merg.append(r1)
            n += 1
            
            # Deal with the last region
            if i == i_max:  # n must > 0
                tab_merg = tab_merg.append(r2)
            else:
                r1 = r2
                ub1 = ub2  # Move the upper boundary to that of r2 for subsequent comparisons

    return tab_merg


if __name__ == '__main__':
    main()
