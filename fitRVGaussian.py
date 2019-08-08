#!/usr/bin/env python3
import argparse, sys, numpy, copy
import matplotlib.pyplot
import datetimelib
import readmultispec
import spectrumlib
import scipy, math
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


class configClass():
	def __init__(self):
		self.range = (0, 1E6)
		self.zoomrange = (0, 1E6)
		self.fitrange = (0, 1E6)
		self.mask = []
		self.ephemFile = None
		self.phaseBins = 10
		self.template = None
		self.normlimits = (-1,-1)
		self.plot = (8.0, 4.9)
		self.title = None
		self.rvfile = None
	
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
				if line[0]=="zoomrange":
					self.zoomrange = (float(line[1]), float(line[2]))
				if line[0]=="restwavelength":
					self.restwavelength = float(line[1])
				if line[0]=="fitrange":
					self.fitrange = (float(line[1]), float(line[2]))
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
				if line[0]=='rvfile':
					self.rvfile = str(line[1])
			except IndexError:
				continue
			except Exception as e:
				print(e)
		fileHandle.close()
		


def sinewave(x, a0, a1, a2):
	period = 1.0
	omega = 2.0 * math.pi / period
	y = a0 + a1 * numpy.sin( omega * (x + a2))
	return y


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Fits a gaussian to an absorption line.')
	parser.add_argument('inputFile', type=str, nargs="+", help='Name of the input file.')
	parser.add_argument('-s', '--save', type=str, help="Save to the plot to a file.")
	parser.add_argument('-p', '--pause', type=float, default=0.5, help="Number of seconds to pause on each plot.")
	parser.add_argument('-i', '--interactive', action="store_true", help="Make each plot interactive (mouse to zoom, etc).")
	parser.add_argument('-x', '--xcor', type=str, help="Optional xcor input data to act as a hint for the starting RV.")
	parser.add_argument('-c', '--config', type=str, help="Configuration file.")
	parser.add_argument('-n', '--noplot', action="store_true", help="Don't plot each fit..")

	arg = parser.parse_args()
	plot=True
	if arg.noplot: plot=False

	blocking = False
	if arg.interactive: blocking = True
	
	config = configClass()
	config.load(arg.config)
	plotWidth = config.plot[0]
	plotHeight = config.plot[1]
	
	if plot:
		spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))
		zoomPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))

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
		print(spectrum)
		spectra.append(spectrum)

	print("Loaded %d spectra"%len(spectra))

	# Trim the spectra
	for s in spectra:
		s.trimWavelengthRange(config.range[0], config.range[1])

	# Attach the phase of the spectra and then sort by phase
	if hasEphem:
		for s in spectra:
			s.phase = ephem.getPhase(s.HJD)
		#spectra = sorted(spectra, key=lambda object: object.phase, reverse = False)
		

	if config.rvfile is not None:
		outfile = open(config.rvfile, 'wt')

	for spectrum in spectra:
		if plot: 
			matplotlib.pyplot.figure(spectrumPlot.number)
			spectrumPlot.canvas.set_window_title(spectrum.name)

			# plot the overall range
			matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'black')
			matplotlib.pyplot.title(spectrum.name)
			matplotlib.pyplot.ylabel(spectrum.fluxUnits)
			matplotlib.pyplot.xlabel(spectrum.wavelengthUnits)
	
			matplotlib.pyplot.figure(zoomPlot.number)
			spectrumPlot.canvas.set_window_title("Zoom ... phase %0.2f"%spectrum.phase)
		zoomedData = copy.deepcopy(spectrum)
		zoomedData.trimWavelengthRange(config.zoomrange[0], config.zoomrange[1])
		if plot: 
			matplotlib.pyplot.step(zoomedData.wavelengths, zoomedData.flux,  color = 'green')
			matplotlib.pyplot.title("%s phase: %0.2f"%(zoomedData.name, zoomedData.phase))
			matplotlib.pyplot.ylabel("Normalised flux")
			matplotlib.pyplot.xlabel("$\mathrm{\AA}$")

		def gaussian(x, a0, a1, a2, a3):
			y = a0 + a1 * numpy.exp(-.5 * ((x-a2)/a3)**2)
			return y

		wavelengthGuess = config.restwavelength + 6 * numpy.sin(2*math.pi*zoomedData.phase)
		print("guess: ", wavelengthGuess)
		fitData = copy.deepcopy(spectrum)
		fitData.trimWavelengthRange(wavelengthGuess - 4.0, wavelengthGuess + 4.0)
		if plot: matplotlib.pyplot.step(fitData.wavelengths, fitData.flux,  color = 'red')

		a0 = 1.0
		a1 = -0.3
		a2 = wavelengthGuess
		a3 = 2.5
		guess = [a0, a1, a2, a3]

		xFit = numpy.arange(min(fitData.wavelengths), max(fitData.wavelengths), 0.1)
		yFit = gaussian(xFit, a0, a1, a2, a3)
		if plot: matplotlib.pyplot.plot(xFit, yFit, linestyle=":")

		xValues = fitData.wavelengths
		yValues = fitData.flux
		yErrors = fitData.fluxErrors
		result, covariance = scipy.optimize.curve_fit(gaussian, xValues, yValues, guess, yErrors)
		errors = numpy.sqrt(numpy.diag(covariance))
		(a0, a1, a2, a3) = result
		yFit = gaussian(xFit, a0, a1, a2, a3)
		
		wavelength = a2
		wavelengthError = errors[2]
		velocity = (wavelength - config.restwavelength)/config.restwavelength * 3E5
		velocityError = 3E5 / config.restwavelength * wavelengthError
		print("%f [%f] A  or %f [%f] km/s"%(wavelength, wavelengthError, velocity, velocityError))
		if config.rvfile is not None:
			outfile.write("%f\t%f\t%f\t\n"%(spectrum.HJD, velocity, velocityError))


		if plot:
			matplotlib.pyplot.plot(xFit, yFit, linestyle="--")
			matplotlib.pyplot.text(6165, 0.85, "Wavelength: %f [%f] $\mathrm{\AA}$"%(wavelength, wavelengthError))
			matplotlib.pyplot.text(6165, 0.80, "Velocity: %f [%f] km/s"%(velocity, velocityError))
			matplotlib.pyplot.text(6165, 0.75, "Width: %f [%f] km/s"%(a3, errors[2]))
		
			matplotlib.pyplot.draw()
			matplotlib.pyplot.show(block=blocking)
			if arg.save is not None:
				print("Writing to file: %s"%arg.save)
				matplotlib.pyplot.savefig(arg.save)
			matplotlib.pyplot.pause(arg.pause)
			spectrumPlot.clf()
			zoomPlot.clf()

	if config.rvfile is not None:
		outfile.close()
