#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import readmultispec
import spectrumlib
from astropy.io import fits


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads FITs spectrum file in IRAF multispec format, plots it and writes out as a JSON file.')
	parser.add_argument('inputFile', type=str, nargs="+", help='Name of the input file.')
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	arg = parser.parse_args()

	plotWidth = 8
	plotHeight = plotWidth/1.62
	spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))

	for index, filename in enumerate(arg.inputFile):
		spectrum = readmultispec.readmultispec(filename)
		hdul = fits.open(filename)
		header = hdul[0]
		specData = header.data
		flux = specData[0][0]
		errors = specData[3][0]
		newSpectrum = spectrumlib.spectrumObject()
		newSpectrum.loadedFromFilename = filename
		newSpectrum.setData(spectrum['wavelen'], flux, errors)
		newSpectrum.extractFITSHeaders(spectrum['header'])
		outputFilename = newSpectrum.name + ".json"
		newSpectrum.writeToJSON(outputFilename)
	

		spectrumPlot.canvas.set_window_title(newSpectrum.name)
		
		matplotlib.pyplot.step(spectrum['wavelen'], spectrum['flux'][0],  color = 'black')

	
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.pause(arg.pause)

		# input("Press enter to continue...")
		spectrumPlot.clf()
		print(" ")

