#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib, configlib
import readmultispec
import spectrumlib, generallib
from astropy.io import fits


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads JSON spectrum file and plots it.')
	parser.add_argument('inputFile', type=str, nargs="+", help='Name of the input file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-c', '--config', type=str, default="plot.cfg", help="Name of the config file. Default is 'plot.cfg'.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	arg = parser.parse_args()

	config = configlib.configClass(debug=False)
	config.load(arg.config)

	blocking = False
	if arg.interactive: blocking = True

	# Set up the matplotlib environment
	generallib.setMatplotlibDefaults()
	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
	matplotlib.rcParams.update(params)

	plotWidth = 7
	plotHeight = plotWidth/1.62
	spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))

	if config.stacked: 
		print("We are going to stack this plot.")
		for index, filename in enumerate(arg.inputFile):
			spectrum = spectrumlib.spectrumObject()
			spectrum.loadFromJSON(filename)
			spectrum.trimWavelengthRange(config.wavelengthRange[0], config.wavelengthRange[1])
			yValues = [f + config.offset*index for f in spectrum.flux]
			matplotlib.pyplot.step(spectrum.wavelengths, yValues,  color = 'black', lw=0.75)
			axes = matplotlib.pyplot.gca()
			axes.set_xlim(config.wavelengthRange[0], config.wavelengthRange[1])
			matplotlib.pyplot.title(config.title)
			matplotlib.pyplot.ylabel(config.ylabel)
			matplotlib.pyplot.xlabel(config.xlabel)
			matplotlib.pyplot.draw()
		
		try: 
			print(config.labels)
			for label in config.labels:
				matplotlib.pyplot.text(label['x'], label['y'], label['text'], size='large')
		except AttributeError:
			print("No labels")
	

		if arg.save is not None:
			print("Writing to file: %s"%arg.save)
			matplotlib.pyplot.savefig(arg.save)
			
		matplotlib.pyplot.show()
		sys.exit()

	
	
	for index, filename in enumerate(arg.inputFile):
		spectrum = spectrumlib.spectrumObject()
		spectrum.loadFromJSON(filename)


		spectrumPlot.canvas.set_window_title(config.title)

		matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'black')
		axes = matplotlib.pyplot.gca()
		axes.set_xlim(config.wavelengthRange[0], config.wavelengthRange[1])
		matplotlib.pyplot.title(config.title)
		matplotlib.pyplot.ylabel(config.ylabel)
		matplotlib.pyplot.xlabel(config.xlabel)

		print(spectrum)
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block=blocking)
		matplotlib.pyplot.pause(arg.pause)
		if arg.save is not None:
			print("Writing to file: %s"%arg.save)
			matplotlib.pyplot.savefig(arg.save)
			

		input("Press enter to continue...")
		spectrumPlot.clf()
		print(" ")

