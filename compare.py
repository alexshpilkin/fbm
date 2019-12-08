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

EPSILON = 1e-6
def close(x, y):
	return abs(x - y) < EPSILON * min(abs(x), abs(y))

def filtervalues(it):
	return filter(lambda line: not line.startswith('#'), it)

def floats(line):
	return [float(x) for x in line.split()]

def readvalues(file1, file2):
	h1, h2 = readheaders(file1), readheaders(file2)
	assert int(h1['Iterations']) == int(h2['Iterations'])
	logn = int(h1['Log of grid size'])
	levels = int(h1['Levels to descend'])

	points, errors = 0, 0
	for line1, line2 in zip(filtervalues(file1), filtervalues(file2)):
		points += 1
		if not all(close(x1, x2)
		           for x1, x2 in zip(floats(line1), floats(line2))):
			errors += 1

	# Wilson score interval for 1 sigma
	return (logn + levels, (errors + 0.5)/(points + 1.0),
	        np.sqrt(errors*(points - errors)/points + 0.25)/(points + 1))

def plotvalues(file1, file2):
	x, y, yerr = readvalues(file1, file2)
	plt.errorbar(x, y, yerr, fmt='.k')

if __name__ == '__main__':
	from sys import argv

	for name1, name2 in zip(argv[1::2], argv[2::2]):
		with open(name1, 'r') as file1, open(name2, 'r') as file2:
			plotvalues(file1, file2)
	plt.gca().xaxis.set_major_locator(MaxNLocator(integer=True))
	plt.xlabel("$L$")
	plt.ylabel("Error probability")
	plt.show()
