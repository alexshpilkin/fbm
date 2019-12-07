#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def readheaders(file):
	headers = {}
	for line in file:
		if line == '#\n': break
		assert line.startswith('# ')
		name, value = line[len('# '):].split(':', maxsplit=1)
		headers[name.strip()] = value.strip()
	return headers

def readerrors(file):
	PREFIX = '# error '

	h = readheaders(file)
	hurst  = float(h['Hurst parameter'])
	logn   = int(h['Log of grid size'])
	levels = int(h['Levels to descend'])

	points = 0
	errors = []
	for line in file:
		if not line.startswith('#'):
			points += 1
			continue
		if not line.startswith(PREFIX):
			continue
		level, *rest = line[len(PREFIX):].split()
		level = levels - int(level)
		errors.append(level)

	bins = np.arange(np.amax(errors) + 1) + logn
	counts = np.bincount(errors)
	return bins, counts/points, np.sqrt(counts)/points

def ploterrors(file, label=None):
	x, y, yerr = readerrors(file)
	plt.errorbar(x, y, yerr, fmt='.', label=label)

if __name__ == '__main__':
	from sys import argv, stdin

	for name in argv[1:]:
		with open(name, 'r') as file:
			ploterrors(file, label=name)
	if len(argv) == 1:
		ploterrors(stdin)
	if len(argv) > 2:
		plt.legend()
	plt.xlim(left=0)
	plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
	plt.xlabel("$\\ell$")
	plt.ylabel("Error probability")
	plt.show()
