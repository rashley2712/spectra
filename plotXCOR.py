#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import readmultispec
import spectrumlib
from astropy.io import fits

class xcorData():
	def __init__(self):
		self.data = []

	def load(self, filename):
		filehandle = open(filename, "rt")
		for line in filehandle:
			if len(line)<2: continue
			if line[0]=="#": continue
			line = line.strip()
			fields = line.split()
			try:
				rv = float(fields[0])
				rv_err = float(fields[1])
				pix = float(fields[2])
				pix_err = float(fields[3])
				integral = float(fields[4])
			except ValueError:
				print("Error parsing line: %s"%line)
			datapoint = {}
			datapoint['rv'] = rv
			datapoint['rv_err'] = rv_err
			datapoint['pix'] = pix
			datapoint['pix_err'] = pix_err
			datapoint['integral'] = integral
			print(datapoint)
			self.data.append(datapoint)

	def getRVs(self):
		rv = [d['rv'] for d in self.data]
		rv_err = [d['rv_err'] for d in self.data]
		return (rv, rv_err)



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads XCOR data from molly and plots it.')
	parser.add_argument('inputFile', type=str, nargs="+", help='Name of the input file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	arg = parser.parse_args()

	blocking = False
	if arg.interactive: blocking = True

	plotWidth = 8
	plotHeight = plotWidth/1.62
	spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))

	for index, filename in enumerate(arg.inputFile):
		xcor = xcorData()
		xcor.load(filename)

		rv, rv_err = xcor.getRVs()
		matplotlib.pyplot.errorbar(numpy.arange(0, 1, 1/len(rv)), rv, marker = '+', color = 'red', yerr = rv_err, lineStyle='none')
		#matplotlib.pyplot.title(spectrum.name)
		#matplotlib.pyplot.ylabel(spectrum.fluxUnits)
		#matplotlib.pyplot.xlabel(spectrum.wavelengthUnits)

		#print(spectrum)
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block=blocking)
		matplotlib.pyplot.pause(arg.pause)
		if arg.save is not None:
			print("Writing to file: %s"%arg.save)
			matplotlib.pyplot.savefig(arg.save)
			

		input("Press enter to continue...")
		spectrumPlot.clf()
		print(" ")

