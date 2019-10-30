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

def readbisects(file):
	PREFIX = '# bisects '

	h = readheaders(file)
	hurst  = float(h['Hurst parameter'])
	logn   = int(h['Log of grid size'])
	levels = int(h['Levels to descend'])

	values = [[] for level in range(levels)]
	for line in file:
		if not line.startswith(PREFIX): continue
		for i, value in enumerate(line[len(PREFIX):].split()):
			values[i].append(int(value))

	level = np.arange(levels) + logn
	length = np.fromiter(map(len, values), int, levels)
	mean = np.fromiter(map(np.mean, values), float, levels)
	std = np.fromiter(map(np.std, values), float, levels)

	return level*hurst, mean, std/np.sqrt(length)

def plotbisects(file, label=None):
	x, y, yerr = readbisects(file)
	plt.errorbar(x, y, yerr, fmt='.', label=label)

if __name__ == '__main__':
	from sys import argv, stdin

	for name in argv[1:]:
		with open(name, 'r') as file:
			plotbisects(file, label=name)
	if len(argv) == 1:
		plotbisects(stdin)
	if len(argv) > 2:
		plt.legend()
	plt.xlim(left=0)
	plt.ylim(bottom=0)
	plt.xlabel("$\ell H$")
	plt.ylabel("Average number of bisections")
	plt.show()
