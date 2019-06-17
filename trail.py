#!/usr/bin/env python3
import argparse, sys, numpy
import matplotlib.pyplot
import datetimelib
import spectrumlib
from astropy.io import fits

class configClass():
	def __init__(self):
		self.lower = 0
		self.upper = 1E6
		self.order = 3
		self.masklower = None
		self.maskupper = None
	
	def load(self, filename="trail.cfg"):
		try:
			fileHandle = open(filename, 'rt')
		except FileNotFoundError:
			print("Could not open config file %s"%filename)
			return
		for line in fileHandle:
			if line[0]=='#': continue
			line = line.strip().split()
			print(line)
			try:
				if line[0]=="wavelength":
					self.lower = float(line[1])
					self.upper = float(line[2])
				if line[0]=="mask":
					self.masklower = float(line[1])
					self.maskupper = float(line[2])
				if line[0]=="order":
					self.order = int(line[1])
			except IndexError:
				continue
			except Exception as e:
				print(e)
		fileHandle.close()
		


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Draws a trailed spectrum from a list of spectra.')
	parser.add_argument('inputFile', type=str, nargs="+", help='Names of the input file.')
	
	arg = parser.parse_args()

	config = configClass()
	config.load()
	
	plotWidth = 8
	plotHeight = plotWidth/1.62
	spectrumPlot = matplotlib.pyplot.figure(figsize=(plotWidth, plotHeight))

	spectra = []

	for index, filename in enumerate(arg.inputFile):
		spectrum = spectrumlib.spectrumObject()
		spectrum.loadFromJSON(filename)
		spectra.append(spectrum)
		print("%d : %s"%(index+1, spectrum.name))
	print("Loaded %d spectra."%len(spectra))

	# Trim to the required wavelength range
	for s in spectra:
		s.trimWavelengthRange(config.lower, config.upper)

	for spectrum in spectra:
		spectrumPlot.canvas.set_window_title(spectrum.name)
		
		matplotlib.pyplot.step(spectrum.wavelengths, spectrum.flux,  color = 'black')
		fit = spectrum.fitPoly(config.order)

		matplotlib.pyplot.plot(spectrum.wavelengths, fit, color='red', linestyle="--")
	
		matplotlib.pyplot.draw()
		matplotlib.pyplot.show(block=False)
		matplotlib.pyplot.pause(0.5)

		key = input("Press enter to continue...")
		if key=='q': sys.exit()
		spectrumPlot.clf()
	
