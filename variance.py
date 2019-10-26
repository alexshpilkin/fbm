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

def readvariances(file):
	PREFIX = '# variance '

	h = readheaders(file)
	hurst  = float(h['Hurst parameter'])
	logn   = int(h['Log of grid size'])
	levels = int(h['Levels to descend'])

	values = [[] for level in range(levels)]
	for line in file:
		if not line.startswith(PREFIX): continue
		level, var = line[len(PREFIX):].split()
		values[levels-int(level)].append(float(var))

	level = np.arange(levels) + logn
	brown = (np.power(2.0, 1-2*hurst) - 0.5) * np.power(2.0, -2*level*hurst)
	length = np.fromiter(map(len, values), int, levels)
	mean = np.fromiter(map(np.mean, values), float, levels)
	std = np.fromiter(map(np.std, values), float, levels)

	return level, mean/brown, std/np.sqrt(length)/brown

def plotvariances(file, fmt='.'):
	x, y, yerr = readvariances(file)
	plt.errorbar(x, y, yerr, fmt=fmt)

if __name__ == '__main__':
	from sys import argv, stdin

	for name in argv[1:]:
		with open(name, 'r') as file:
			plotvariances(file)
	if len(argv) == 1:
		plotvariances(stdin)
	if len(argv) > 2:
		plt.legend(argv[1:])
	plt.show()
