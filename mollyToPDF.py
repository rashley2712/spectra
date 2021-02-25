#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import spectrumlib, generallib, configlib
import trm.molly
from matplotlib.backends.backend_pdf import PdfPages


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Loads Molly (.mol) file of spectra and plots them all as a PDF file. It needs a configuration file call plot.cfg.')
	parser.add_argument('inputFile', type=str, help='Name of the molly file.')
	parser.add_argument('outputFile', type=str, help="The filename of the pdf file.")
	parser.add_argument('-c', '--config', type=str, default='plot.cfg', help="Config file. Default is 'plot.cfg'")
	parser.add_argument('-n', '--number', type=int, default=0, help="Number of the spectrum in the molly file (starting with 1). Default is 0 = 'all'.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Show an interactive preview after each page (mouse to zoom, etc).")
	parser.add_argument('-p', '--plotsPerPage', type=int, default=3, help="Number of plots per page. Default is 3.")
	arg = parser.parse_args()

	# Set up the configuration
	config = configlib.configClass(debug=False)
	config.load(filename=arg.config)
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
		if spectrum.fluxUnits=="MILLIJANSKYS": spectrum.fluxUnits = "mJy"
		# spectrum.fluxUnits = "relative counts"
		
		print("Parsed headers of %s for HJD: %f"%(targetName, spectrum.HJD))
		spectrum.name = "%s-%f"%(spectrum.objectName, spectrum.HJD)
		spectra.append(spectrum)
		print(spectrum)
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
	
	
	if arg.number>0:
		print("Plotting spectrum number: %s"%arg.number)
		if arg.number>numSpectra:
			print("Spectrum out of range. There are %d spectra in the molly file %s."%(numSpectra, arg.inputFile))
			sys.exit()
		spectra = [spectra[arg.number-1]]
		plotsPerPage = 1
	else:
		print("Plotting all spectra %d in the file."%numSpectra)
		plotsPerPage = arg.plotsPerPage

	with PdfPages(arg.outputFile) as pdf:
		pageNumber = 1
		for index, spectrum in enumerate(spectra):
			colIndex = index % plotsPerPage + 1
			if colIndex == 1:
				spectrumPlot = matplotlib.pyplot.figure(figsize=(8, 11))
			print("page number: %d  plot number: %d"%(pageNumber, colIndex))
			matplotlib.pyplot.subplot(plotsPerPage, 1, colIndex )
			spectrumPlot.tight_layout(pad=2.0) 
			if hasattr(config, 'title'):
				spectrumPlot.canvas.set_window_title(config.title)
				if config.title != "None":  matplotlib.pyplot.title(config.title)
			else:	
				spectrumPlot.canvas.set_window_title(spectrum.name)
				matplotlib.pyplot.title(spectrum.name)

			matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'black', lw=0.5)
			axes = matplotlib.pyplot.gca()
			axes.set_xlim(config.wavelengthRange[0], config.wavelengthRange[1])
			matplotlib.pyplot.ylabel(config.ylabel + " (" + spectrum.fluxUnits + ")")
			if (colIndex == plotsPerPage): 
				matplotlib.pyplot.draw()
				matplotlib.pyplot.xlabel(config.xlabel)
				pdf.savefig()
				pageNumber+=1
				if arg.interactive: matplotlib.pyplot.show(block=True)
		if colIndex!=plotsPerPage:
			pdf.savefig()
	matplotlib.pyplot.close()

	
