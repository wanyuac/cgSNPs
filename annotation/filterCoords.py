#!/bin/env python

"""
This script is a Python3-compatible version of script filterCoords.py in the RedDog package
(https://github.com/katholt/RedDog). It filters Mummer inexact repeat-match coord output to coordinate
table for parseSNPTable.py.

Example command: python filterCoords.py -i seq_seq.coords -o seq_seq_filtered.coords -I 90

RedDog:
	Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
	All rights reserved. (see README.txt for more details)
	
Created: 10 Sep 2014
Changes:
         24 Apr 2020 - Update for Python3 syntax - Yu Wan <wanyuac@126.com>
"""

import os, sys, glob
import subprocess
import string
import re
import operator
from optparse import OptionParser


def main():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	parser.add_option("-i", "--input", action="store", dest="input", help="Input: coords file from Mummer numcer (default none)", default="")
	parser.add_option("-o", "--output", action="store", dest="output", help="Output: filtered coords file for parseSNPTable (default none)", default="")
	parser.add_option("-I", "--Identity", action="store", dest="Identity", help="Identity: level of Identity to filter out (default >=85.0)", default=85.0)
	
	return parser.parse_args()


if __name__ == "__main__":

	(options, args) = main()
	
	def filterCoords(coords, input, Identity):

		in_file = open(input)
		count = 0
		for line in in_file:
			if count <= 5:
				count += 1
			else:
				data = line.split("|")
				if float(data[3]) >= Identity:
					new_coords = data[0].split()
					start = int(new_coords[0])
					stop = int(new_coords[1])
					coords.append([str(start),str(stop)])		
		in_file.close()
		return(coords)

	def coordsOut(coords, output):
		coords_out = ""
		count = 0
		for coord in coords:
			coords_out += (coord[0] + "," + coord[1] +"\n")
		output_file = open(output, "w")
		output_file.write(coords_out)
		output_file.close()
		return
		
	### MAIN PROCESS
			
	# set up variables
	if options.input == '':
		print('\nNo input repeats file provided (-i)')
		sys.exit()
	else:
		input = options.input
	if options.output == '':
		print('\nNo output coords file provided (-o)')
		sys.exit()
	else:
		output = options.output
	Identity = float(options.Identity)
	coords = []
	
	coords = filterCoords(coords, input, Identity)
	coordsOut(coords, output)
