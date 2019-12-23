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

	# Wilson score interval for 1 sigma
	return (epsilon, (errors + 0.5)/(points + 1.0),
	        np.sqrt(errors*(points - errors)/points + 0.25)/(points + 1))

def plotepsilon(file, fmt='.k'):
	x, y, yerr = readepsilon(file)
	plt.errorbar(x, y, yerr, fmt=fmt)

if __name__ == '__main__':
	from sys import argv, stdin

	fmt = '.k'
	for name in argv[1:]:
		if name.startswith('-f'):
			fmt = None
			continue
		if fmt is None:
			fmt = name
			continue
		with open(name, 'r') as file:
			plotepsilon(file, fmt=fmt)
	if len(argv) == 1:
		plotepsilon(stdin, fmt=fmt)
	plt.plot(plt.xlim(), np.array(plt.xlim()) * 3, '--', color='gray', linewidth=.5)
	plt.plot(plt.xlim(), np.array(plt.xlim()) * 10, '--', color='gray', linewidth=.5)
	plt.xscale('log')
	plt.xlabel("$\epsilon'$")
	plt.yscale('log')
	plt.ylabel("Error probability")
	plt.show()
