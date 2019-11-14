#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import spectrumlib, generallib, configlib
import trm.molly


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads Molly (.mol) file of spectra and plots them. It needs a configuration file.')
	parser.add_argument('inputFile', type=str, help='Name of the molly file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	arg = parser.parse_args()

	# Set up the configuration
	config = configlib.configClass(debug=False)
	config.load()
	config.show()

	blocking = False
	if arg.interactive: blocking = True

	# Set up the matplotlib environment
	generallib.setMatplotlibDefaults()
	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
	matplotlib.rcParams.update(params)

	
	mollyFilename = arg.inputFile
	mollyFile = trm.molly.rmolly(mollyFilename)
    
	spectra = []
		
	for index, r in enumerate(mollyFile):
		wavelengths = []
		flux = []
		fluxErrors = []	
		for f, fe, w in zip(r.f, r.fe, r.wave):
			# print(w, f, fe)
			wavelengths.append(w)
			flux.append(f)
			fluxErrors.append(fe)
 
		head = r.head
		spectrum = spectrumlib.spectrumObject()
		npoints = spectrum.setData(wavelengths, flux, fluxErrors)
		targetName = spectrum.parseHeaderInfo(head)
		spectrum.wavelengthUnits = "A"
		spectrum.fluxLabel = r.label
		spectrum.fluxUnits = r.units
		# spectrum.fluxUnits = "relative counts"
		
		print("Parsed headers of %s for HJD: %f"%(targetName, spectrum.HJD))
		spectrum.name = "%s-%f"%(spectrum.objectName, spectrum.HJD)
		spectra.append(spectrum)
		
	numSpectra = len(spectra)

	print("%d spectra loaded."%numSpectra)

	if hasattr(config, 'wavelengthRange'):
		print("Trimming wavelength range to (%f, %f)."%(config.wavelengthRange[0], config.wavelengthRange[1]))
		for s in spectra:
			s.trimWavelengthRange(config.wavelengthRange[0], config.wavelengthRange[1])

	if hasattr(config, 'plotDimensions'):
		plotWidth, plotHeight = config.plotDimensions
	else:
		plotWidth = 8
		plotHeight = plotWidth/1.62
	
	spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))

	print("Plotting spectrum number: %s", config.number)
	if config.number>0:
		spectra = [spectra[config.number]]

	for spectrum in spectra:
		
		if hasattr(config, 'title'):
			spectrumPlot.canvas.set_window_title(config.title)
			if config.title != "None":  matplotlib.pyplot.title(config.title)
		else:	
			spectrumPlot.canvas.set_window_title(spectrum.name)
			matplotlib.pyplot.title(spectrum.name)

		matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'black', lw=0.5)
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
	print(" ")

