#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import readmultispec
import spectrumlib
from astropy.io import fits


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads JSON spectrum file and plots it.')
	parser.add_argument('inputFile', type=str, nargs="+", help='Name of the input file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	arg = parser.parse_args()

	plotWidth = 8
	plotHeight = plotWidth/1.62
	spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))

	for index, filename in enumerate(arg.inputFile):
		spectrum = spectrumlib.spectrumObject()
		spectrum.loadFromJSON(filename)

		spectrumPlot.canvas.set_window_title(spectrum.name)

		matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'black')
		matplotlib.pyplot.title(spectrum.name)
		matplotlib.pyplot.ylabel(spectrum.fluxUnits)
		matplotlib.pyplot.xlabel(spectrum.wavelengthUnits)

		print(spectrum)
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.pause(arg.pause)
		if arg.save is not None:
			print("Writing to file: %s"%arg.save)
			matplotlib.pyplot.savefig(arg.save)
			

		input("Press enter to continue...")
		spectrumPlot.clf()
		print(" ")

