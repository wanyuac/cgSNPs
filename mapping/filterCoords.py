#!/bin/env python

"""
This script is a Python3-compatible version of script filterCoords.py in the RedDog package
(https://github.com/katholt/RedDog). It filters Mummer inexact repeat-match coord output to coordinate
table for parseSNPTable.py.

Example command: python filterCoords.py -i seq_seq.coords -o seq_seq_filtered.coords -I 90

RedDog:
	Copyright (c) 2015, David Edwards, Bernie Pope, Kat Holt
	All rights reserved. (see README.txt for more details)

Python v2 and v3 (recommended) compatible.

Created: 10 Sep 2014
Changes:
		24 Apr 2020 - Update for Python3 syntax - Yu Wan <wanyuac@126.com>
		28 Apr 2020 - Fix a critical issue that the script failed to ignore self-hits - Yu Wan
		1 May 2020 - Enable output to stdout and rewrote the script - Yu Wan
Latest update: 9 Jan 2022
"""
from __future__ import print_function  # To enable Python 2 to interpret the command print(line, end = "").
from argparse import ArgumentParser


def parse_arguments():
	parser = ArgumentParser(description = "Filter coordinates from program show-coords by nucleotide identity")
	parser.add_argument("-i", "--input", type = str, dest = "input", required = True,\
						help = "Input: coords file from Mummer numcer (default none)")
	parser.add_argument("-o", "--output", type = str, dest = "output", required = False, default = "",\
						help = "Output: filtered coords file for parseSNPTable (default stdout)")
	parser.add_argument("-I", "--Identity", type = float, dest = "Identity", required = False, default = 90.0,\
						help = "Identity: level of Identity to filter out (default >=90.0)")
	
	return parser.parse_args()


def main():
	arguments = parse_arguments()
	coordsOut(filterCoords(arguments.input, arguments.Identity), arguments.output)
		
	return


def filterCoords(input, Identity):
	coords = []
	count = 0
	
	in_file = open(input, "r")
	for line in in_file:
		if count <= 5:  # Skip the first six lines that do not provide hit information
			count += 1
		else:  # Start processing hits from the seventh line of the input file.
			data = line.split("|")
			if float(data[3]) >= Identity:
				coords_ref = data[0].split()  # Drop all white-spaces
				coords_query = data[1].split()
				"""
				S1, E1: with respect to the reference sequence; S2, E2: with respect to the query
				sequence. See http://mummer.sourceforge.net/manual.
				"""
				start_r = int(coords_ref[0])
				stop_r = int(coords_ref[1])
				start_q = int(coords_query[0])
				stop_q = int(coords_query[1])
				
				# Skip all self-hits
				if start_r != start_q or stop_r != stop_q:
					coords.append([str(start_r),str(stop_r)])
	in_file.close()
	
	return(coords)


def coordsOut(coords, output):
	coords_out = ""
	count = 0
	for coord in coords:
		coords_out += (coord[0] + "," + coord[1] +"\n")
	
	if output == "":
		for line in coords_out:
			print(line, end = "")  # Do not append an extra newline character to every line. (A feature introduced by Python 3)
	else:
		with open(output, "w") as output_file:
			output_file.write(coords_out)
		
	return


if __name__ == "__main__":
	main()
