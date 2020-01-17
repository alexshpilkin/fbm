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

def readbridges(file):
	PREFIX = '# bridge '

	h = readheaders(file)
	logn   = int(h['Log of grid size'])
	levels = int(h['Levels to descend'])

	levelhists = []
	timehists  = []
	poshists   = []
	levelhist = []
	timehist  = []
	poshist   = []
	for line in file:
		if not line.startswith('#'):
			levelhists.append(levelhist)
			levelhist = []
			timehists.append(timehist)
			timehist = []
			poshists.append(poshist)
			poshist = []
		if not line.startswith(PREFIX):
			continue
		level, ltime, lpos, rtime, rpos = line[len(PREFIX):].split()
		levelhist.append(logn + levels - int(level))
		timehist.append([float(ltime), float(rtime)])
		poshist.append([float(lpos), float(rpos)])

	return levelhists, timehists, poshists

def plotbridges(file, index, label=None):
	levelhists, timehists, poshists = readbridges(file)
	levelhist, timehist, poshist = levelhists[index], timehists[index], poshists[index]
	line, = plt.plot(timehist[0], poshist[0], '-',
	                 label=label, linewidth=10.0/levelhist[0])
	for level, time, pos in zip(levelhist[1:], timehist[1:], poshist[1:]):
		plt.plot(time, pos, '-',
		         color=line.get_color(), linewidth=10.0/level)

if __name__ == '__main__':
	from sys import argv, stdin
	for name, index in zip(argv[1::2], argv[2::2]):
		with open(name, 'r') as file:
			plotbridges(file, int(index),
			            label='{}[{}]'.format(name, index))
	if len(argv) == 2:
		plotbridges(stdin, int(argv[1]))
	if len(argv) > 3:
		plt.legend()
	plt.xlabel("Time")
	plt.ylabel("Position")
	plt.show()
