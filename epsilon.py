#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def readheaders(file):
	headers = {}
	for line in file:
		if line == '#\n': break
		assert line.startswith('# ')
		name, value = line[len('# '):].split(':', maxsplit=1)
		headers[name.strip()] = value.strip()
	return headers

def readepsilon(file):
	PREFIX = '# error '

	h = readheaders(file)
	hurst = float(h['Hurst parameter'])
	epsilon = float(h['Error tolerance'])

	points = 0
	errors = 0
	for line in file:
		if not line.startswith('#'):
			points += 1
			continue
		if not line.startswith(PREFIX):
			continue
		errors += 1

	return epsilon, errors/points, np.sqrt(errors)/points

def plotepsilon(file):
	x, y, yerr = readepsilon(file)
	plt.errorbar(x, y, yerr, fmt='.k')

if __name__ == '__main__':
	from sys import argv, stdin

	for name in argv[1:]:
		with open(name, 'r') as file:
			plotepsilon(file)
	if len(argv) == 1:
		plotepsilon(stdin)
	plt.xscale('log')
	plt.xlabel("$\epsilon'$")
	plt.yscale('log')
	plt.ylabel("Error probability")
	plt.show()
