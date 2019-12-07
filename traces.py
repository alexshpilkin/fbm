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

def readtraces(file):
	PREFIX = '# variance '

	h = readheaders(file)
	hurst  = float(h['Hurst parameter'])
	logn   = int(h['Log of grid size'])
	levels = int(h['Levels to descend'])

	levelhists = []
	levelhist  = []
	varhists = []
	varhist  = []
	for line in file:
		if not line.startswith('#'):
			levelhists.append(levelhist)
			levelhist = []
			varhists.append(varhist)
			varhist = []
		if not line.startswith(PREFIX):
			continue
		level, var = line[len(PREFIX):].split()
		level = logn + levels - int(level)
		brown = (np.power(2.0, 1-2*hurst) - 0.5) * np.power(2.0, -2*int(level)*hurst)
		levelhist.append(level)
		varhist.append(float(var) / brown)

	return 2**logn, levelhists, varhists

def plottraces(file, index, label=None):
	n, levelhists, varhists = readtraces(file)
	plt.subplot(2, 1, 1)
	plt.plot(np.arange(n, len(levelhists[index])+n), levelhists[index],
	         label=label)
	plt.subplot(2, 1, 2)
	plt.plot(np.arange(n, len(varhists[index])+n), varhists[index])

if __name__ == '__main__':
	from sys import argv, stdin
	for name, index in zip(argv[1::2], argv[2::2]):
		with open(name, 'r') as file:
			plottraces(file, int(index),
			           label='{}[{}]'.format(name, index))
	if len(argv) == 2:
		plottraces(stdin, int(argv[1]))
	plt.subplot(2, 1, 1)
	plt.ylim(bottom=0)
	plt.ylabel("$\ell$")
	if len(argv) > 3:
		plt.legend()
	plt.subplot(2, 1, 2)
	plt.ylabel("$\\sigma^2 / 2^{-2\\ell H}(2^{1-2H}-1/2)$")
	plt.xlabel("Number of points sampled")
	plt.show()
