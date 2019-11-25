#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from mpmath import fp
from scipy.special import erfi

def hyp2f2(a, b, c, d, zs):
	return np.vectorize(lambda z: fp.hyp2f2(a, b, c, d, z))(zs)

def eye(z):
	z = np.asarray(z, float)
	return (z**4/6 * hyp2f2(1.0, 1.0, 2.5, 3, z*z/2) +
	        np.pi * (1 - z*z) * erfi(z / np.sqrt(2)) - 3 * z*z +
	        np.sqrt(2*np.pi) * np.exp(z*z/2) * z + 2.0)

def f1(y):
	y = np.asarray(y, float)
	return eye(y) + y*y * (np.log(2 * y*y) + np.euler_gamma) - 2.0

def maxprob(hurst, y):
	y = np.asarray(y, float)
	return (np.power(y, 1.0/hurst - 2) / (np.sqrt(2*np.pi) * hurst) *
	        np.exp(-y*y/2 + (hurst-0.5) * (f1(y) - 2*(np.euler_gamma + np.log(2)))))

def readheaders(file):
	headers = {}
	for line in file:
		if line == '#\n': break
		assert line.startswith('# ')
		name, value = line[len('# '):].split(':', maxsplit=1)
		headers[name.strip()] = value.strip()
	return headers

def readpositions(file):
	h = readheaders(file)
	hurst = float(h['Hurst parameter'])
	TIME = 1.0

	data = []
	for line in file:
		if line.startswith('#'): continue
		data.append(float(line) / np.sqrt(2) / TIME**hurst)

	nhist, bins = np.histogram(data, 50, density=True)
	uhist, _    = np.histogram(data, 50)
	y = (bins[:-1] + bins[1:])/2
	pth = maxprob(hurst, y)

	yy = np.arange(0.005, 8, 0.005)
	pth /= np.trapz(maxprob(hurst, yy), yy)

	return y[1:], nhist[1:], (nhist / np.sqrt(uhist))[1:], pth[1:]

def plotpositions(file, label=None):
	y, p, perr, pth = readpositions(file)
	line, = plt.plot(y, pth, label=label)
	plt.errorbar(y, p, perr, fmt='.', color=line.get_color())

if __name__ == '__main__':
	from sys import argv, stdin

	for name in argv[1:]:
		with open(name, 'r') as file:
			plotpositions(file, label=name)
	if len(argv) == 1:
		plotpositions(stdin)
	if len(argv) > 2:
		plt.legend()
	plt.xlim(left=0)
	plt.ylim(bottom=0)
	plt.xlabel("$y$")
	plt.ylabel("$P$")
	plt.show()
