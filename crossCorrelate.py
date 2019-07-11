#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import spectrumlib
from astropy.io import fits


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Cross correlates a spectrum(or spectra) with a template.')
	parser.add_argument('inputFile', type=str, help='Name of the input file.')
	parser.add_argument('templateFile', type=str, help='Name of the template file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	arg = parser.parse_args()
	
	blocking = False
	if arg.interactive: blocking = True
    
	plotWidth = 8
	plotHeight = plotWidth/1.62
	
	templateSpectrum = spectrumlib.spectrumObject()
	templateSpectrum.loadFromJSON(arg.templateFile)
	print("Template spectrum:")
	print(templateSpectrum)
	print()
	templatePlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	templatePlot.canvas.set_window_title("template")
	matplotlib.pyplot.step(templateSpectrum.wavelengths, templateSpectrum.flux,  color = 'green')
	matplotlib.pyplot.title(templateSpectrum.name)
	matplotlib.pyplot.ylabel(templateSpectrum.fluxUnits)
	matplotlib.pyplot.xlabel(templateSpectrum.wavelengthUnits)

	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block=False)

	spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	spectrum = spectrumlib.spectrumObject()
	spectrum.loadFromJSON(arg.inputFile)
	print("Target spectrum:")
	print(spectrum)
	print()
	matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'red')
	matplotlib.pyplot.title(spectrum.name)
	matplotlib.pyplot.ylabel(spectrum.fluxUnits)
	matplotlib.pyplot.xlabel(spectrum.wavelengthUnits)

	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block=False)
	matplotlib.pyplot.pause(arg.pause)

	target = spectrum.flux
	reference = templateSpectrum.flux
	cross = numpy.correlate(target, reference)
	print(cross)
	

	from scipy import signal
	corr = signal.correlate(target, reference, mode='full')
	print(corr, len(corr))
	for c in corr:
		print(c)
	max = numpy.max(corr)
	index = numpy.argmax(corr)
	print("Max: %f at %f"%(max, index))

	correlationPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
	matplotlib.pyplot.plot(range(0, len(corr)), corr, color = 'black')
	matplotlib.pyplot.title("Cross correlation...")
	#matplotlib.pyplot.ylabel(spectrum.fluxUnits)
	#matplotlib.pyplot.xlabel(spectrum.wavelengthUnits)

	matplotlib.pyplot.draw()
	matplotlib.pyplot.show(block=True)
	matplotlib.pyplot.pause(arg.pause)
	sys.exit()

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
		matplotlib.pyplot.show(block=blocking)
		matplotlib.pyplot.pause(arg.pause)
		if arg.save is not None:
			print("Writing to file: %s"%arg.save)
			matplotlib.pyplot.savefig(arg.save)
			

		input("Press enter to continue...")
		spectrumPlot.clf()
		print(" ")

	input("Press enter to continue...")
		