#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import readmultispec
import spectrumlib
from astropy.io import fits


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Reads a file of HJD and ephemeris and computes a RV for each line.')
	parser.add_argument('inputFile', type=str, help='Name of the input file.')
	parser.add_argument('-e', '--ephemeris', type=str, default='ephemeris.dat', help='Ephemeris file (default: ephemeris.dat).')
	parser.add_argument('-o', '--outputFile', type=str, default='rvs.dat', help='Output file.')
	
	arg = parser.parse_args()
	ephem = datetimelib.ephemerisObject()
	ephem.loadFromFile(arg.ephemeris)
	print(ephem)
	

	inputFile = open(arg.inputFile, "rt")	
	outputFile = open(arg.outputFile, "wt")

	for line in inputFile:
		if line[0]=="#": continue
		if len(line)<3: continue
		line = line.strip()
		fields = line.split()
		print(fields)
		id = int(fields[0])
		#run = int(fields[1])
		HJD = float(fields[1])
		print("HJD: %f   velocity: %f"%(HJD, ephem.getRV(HJD)))
		outputFile.write("%f\n"%(ephem.getRV(HJD)))
		

	inputFile.close()

	outputFile.close()

	
