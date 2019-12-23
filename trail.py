#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import spectrumlib
import generallib
from astropy.io import fits

class configClass():
	def __init__(self):
		self.range = (0, 1E6)
		self.mask = []
		self.ephemFile = None
		self.phaseBins = 10
		self.template = None
		self.normlimits = (-1,-1)
		self.plot = (8.0, 4.9)
		self.title = None
	
	def load(self, filename="trail.cfg"):
		try:
			fileHandle = open(filename, 'rt')
		except FileNotFoundError:
			print("Could not open config file %s"%filename)
			return
		for line in fileHandle:
			if line[0]=='#': continue
			original = line
			line = line.strip().split()
			print(line)
			try:
				if line[0]=="range":
					self.range = (float(line[1]), float(line[2]))
				if line[0]=="mask":
					self.mask.append((float(line[1]), float(line[2])))
				if line[0]=="order":
					self.order = int(line[1])
				if line[0]=="ephem":
					self.ephemFile = str(line[1])
				if line[0]=="title":
					self.title = str(original[6:-1])
					print("Title: ", self.title)
				if line[0]=="phasebins":
					self.phaseBins = int(line[1])
				if line[0]=="template":
					self.template = str(line[1])
					self.normlimits = (float(line[2]), float(line[3]))
				if line[0]=="plot":
					self.plot = (float(line[1]), float(line[2]))
			except IndexError:
				continue
			except Exception as e:
				print(e)
		fileHandle.close()
		


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Draws a trailed spectrum from a list of spectra.')
	parser.add_argument('inputFile', type=str, nargs="+", help='Names of the input file.')
	parser.add_argument('--noplot', action="store_true", help="Suppress plotting of the input spectra.")	
	parser.add_argument('-c', '--config', type=str, help="Configuration file for the trail plots. Default is 'trail.cfg'.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-b', '--boost', action="store_true", help="Boost the trail contrast.")
	
	# Set up the matplotlib environment
	generallib.setMatplotlibDefaults()
	params = {	'axes.labelsize': 'large',
				'xtick.labelsize': 'large',
				'ytick.labelsize': 'large',
			}
	matplotlib.rcParams.update(params)
	
	arg = parser.parse_args()
	plot = True
	if arg.noplot:
		plot = False

	config = configClass()
	config.load(arg.config)
	plotWidth = config.plot[0]
	plotHeight = config.plot[1]

	hasEphem = False
	if config.ephemFile is not None:
		ephem = datetimelib.ephemerisObject()
		ephem.loadFromFile(config.ephemFile)
		print(ephem)
		hasEphem = True
	
	spectra = []

	for index, filename in enumerate(arg.inputFile):
		spectrum = spectrumlib.spectrumObject()
		spectrum.loadFromJSON(filename)
		spectra.append(spectrum)
		print("%d : %s"%(index+1, spectrum.name))
	print("Loaded %d spectra."%len(spectra))

	# Attach the phase of the spectra and then sort by phase
	if hasEphem:
		for s in spectra:
			s.phase = ephem.getPhase(s.HJD)
		spectra = sorted(spectra, key=lambda object: object.phase, reverse = False)
		


	if config.template is not None:
		# Perform the normalisation across all spectra
		referenceSpectrum = None
		for s in spectra:
			if s.name == config.template: referenceSpectrum = s

		if referenceSpectrum is None:
			referenceSpectrum = spectra[0]
		print("Reference spectrum: \n%s"%referenceSpectrum)
		normalConstant = referenceSpectrum.integrate(config.normlimits)
		print("Normalisation constant: %f"%normalConstant)
	
		for s in spectra:
			normalVal = s.integrate(config.normlimits)
			print("Normalisation value:", normalVal, normalConstant)
			#s.divide(normalVal/normalConstant)
			s.divide(normalConstant)
	

#	else:
		#continue
		# Normalise original spectrum by dividing by the fit 
		#for s in spectra:
		#	fit = s.fitPoly(config.order, config.mask)
		#	s.divideArray(fit)		

	# Trim to the required wavelength range
	for s in spectra:
		print(s.wavelengthRange)
		print(len(s.wavelengths), len(s.flux), len(s.fluxErrors))
		if len(s.fluxErrors)==0: s.fluxErrors = numpy.ones(len(s.wavelengths))
		s.trimWavelengthRange(config.range[0], config.range[1])

	# Resample down to lowest resolution
	minElements = 1E6
	for index, s in enumerate(spectra):
		elements = len(s.flux)
		if elements < minElements:
			minElements = elements
			shortest = index
	print("Shortest spectrum is number %d with %d elements."%(shortest, minElements))
	sampleWavelengths = spectra[shortest].wavelengths
	for s in spectra:
		s.resample(sampleWavelengths)
	

	# Remove any negative values in the spectra
	for s in spectra:
		spectrum.removeNegatives()
	
	if plot:
		spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		normalisedPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		for spectrum in spectra:
			matplotlib.pyplot.figure(spectrumPlot.number)
			spectrumPlot.canvas.set_window_title(spectrum.name)
			matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'black')
	
		 
		#	if config.template is None: matplotlib.pyplot.plot(spectrum.wavelengths, fit, color='red', linestyle="--")
			axes = matplotlib.pyplot.gca()

			for m in config.mask:
				width = m[1] - m[0]
				height = axes.get_ylim()[1] - axes.get_ylim()[0]
				rect = matplotlib.pyplot.Rectangle((m[0],axes.get_ylim()[0]),width, height, linewidth=1, edgecolor=None,fill=False, facecolor='blue', hatch='//')
				axes.add_patch(rect) 

			
			matplotlib.pyplot.draw()
		
			matplotlib.pyplot.figure(normalisedPlot.number)
			normalisedPlot.canvas.set_window_title("Normalised: " + spectrum.name)
			matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'black')
			
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show(block=False)
			matplotlib.pyplot.pause(arg.pause)

			key = ""
			#key = input("Press enter to continue...")

			if key=='q': sys.exit()
			spectrumPlot.clf()
			normalisedPlot.clf()

	


	ySize = len(spectra)
	xSizeArray = [len(s.wavelengths) for s in spectra]
	xSize = max(xSizeArray)
	trailArray = []
	if not hasEphem:
		for index, s in enumerate(spectra):
			trailArray.append(s.flux)
			print(len(s.flux), s.name, s.loadedFromFilename)
		trailBitmap = numpy.copy(trailArray)
	else:
		numPhaseBins = config.phaseBins
		ySize = numPhaseBins * 2
		xSizeArray = [len(s.wavelengths) for s in spectra]
		xSize = max(xSizeArray)
		trailArray = []
		phasedArray = []
		for index in range(ySize):
			blankSpectrum = numpy.zeros(xSize)
			trailArray.append(blankSpectrum)
			phasedArray.append(blankSpectrum)
		
		for bin in range(ySize):
			phaseLower = 2.0/ySize * bin
			phaseUpper = 2.0/ySize * (bin+1)
			spectraToAdd = []
			for s in spectra:
				if s.phase<phaseUpper and s.phase>=phaseLower: 
					spectraToAdd.append(s)
				if (s.phase+1)<phaseUpper and (s.phase+1)>=phaseLower: 
					spectraToAdd.append(s)
			totalSpectrum = numpy.zeros(xSize)
			for s in spectraToAdd:
				totalSpectrum+=s.flux
			averageSpectrum = totalSpectrum / len(spectraToAdd)
			phasedArray[bin] = averageSpectrum
			print("Bin #%d  phase range %f to %f contains %d spectra to average together."%(bin, phaseLower, phaseUpper, len(spectraToAdd)))
			

	
		trailBitmap = numpy.array(numpy.copy(phasedArray))

	startWavelength = min(s.wavelengths)
	endWavelength = max(s.wavelengths)
	
	if hasEphem: yMax = 2
	else: yMax = len(spectra)
	
	if arg.boost: trailBitmap = generallib.percentiles(trailBitmap, 20, 99)
	
	trailPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))	
	
	#matplotlib.pyplot.imshow(trailBitmap, extent=[startWavelength, endWavelength, yMax, 0])
	wavelength2velocity = lambda w: (w - 6562.8) / 6562.8 * 3E5 
	
	#matplotlib.pyplot.imshow(trailBitmap, extent=[wavelength2velocity(startWavelength), wavelength2velocity(endWavelength), yMax, 0])
	matplotlib.pyplot.imshow(trailBitmap, extent=[startWavelength, endWavelength, yMax, 0])
	#matplotlib.pyplot.imshow(trailBitmap)
	axes = matplotlib.pyplot.gca()
	axes.invert_yaxis()
	matplotlib.pyplot.axis('auto')
	axes.figure.set_size_inches(plotWidth, plotHeight)
	
	#newPositions = [0, 10, 20, 200]
	#axes.set_xticks(newPositions)

	# this is code to add a velocity scale to to the top of the x-axis
	#ax2 = axes.twiny()
	#wavelength2velocity = lambda w: (w - 6562.8) / 6562.8 * 3E5 
	#newlabel = [273.15,290,310,330,350,373.15] # labels of the xticklabels: the position in the new x-axis
	#newLabels = [-2500, 0, 2500]
	#newPositions = [wavelength2velocity(w) for w in s.wavelengths]
	#print(newPositions)
	#k2degc = lambda t: t-273.15 # convert function: from Kelvin to Degree Celsius
	#newpos   = [k2degc(t) for t in newlabel]   # position of the xticklabels in the old x-axis
	#ax2.set_xticks(newPositions)
	#ax2.set_xticklabels(newLabels)
	#ax2.set_xlabel('Velocity [km/s]')
	# ax2.set_xlim(axes.get_xlim())

	if config.title is not None:
		matplotlib.pyplot.title(config.title)
	
	if hasEphem: matplotlib.pyplot.ylabel('Phase')
	else: matplotlib.pyplot.ylabel('Spectrum number')
	matplotlib.pyplot.xlabel('Wavelength (\AA)')
	
	spectScale = (endWavelength - startWavelength) / xSize
	spectrum2velocity = lambda n: wavelength2velocity((n*spectScale)+startWavelength)
	velocityToWavelength = lambda v: v / 3E5 * 6562.8 + 6562.8
	topAxes = axes.twiny()
	velocities = [spectrum2velocity(n) for n in range(xSize)]
	#velocities = [wavelength2velocity(w) for w in spectrum.wavelengths]
	#print(velocities)
	bottomLabels = axes.get_xticklabels()
	print(bottomLabels)
	velocityLabels = [-1500, -1000, -500,  0, 500, 1000, 1500]
	xTicks = [velocityToWavelength(v) for v in velocityLabels]
	print(velocityLabels)
	print(xTicks)
	topAxes.set_xticks(xTicks)	
	topAxes.set_xticklabels(velocityLabels)	
	print(topAxes.get_xlim())
	topAxes.set_xlim(axes.get_xlim())
	topAxes.set_xlabel('Velocity ($\mathrm{km}\,\mathrm{s}^{-1}$)')
	
	matplotlib.pyplot.draw()
	if arg.save is not None:
			print("Writing to file: %s"%arg.save)
			matplotlib.pyplot.savefig(arg.save)
	matplotlib.pyplot.show(block=True)

