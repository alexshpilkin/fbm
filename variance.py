#!/usr/bin/env python3

from sys import argv, exit, stderr, stdin
import numpy as np
import matplotlib.pyplot as plt

if len(argv) != 4:
	print("usage: {} HURST LOGN LEVELS".format(argv[0]), file=stderr)
	exit(64) # EX_USAGE
hurst  = float(argv[1])
logn   = int(argv[2])
levels = int(argv[3])

values = [[] for level in range(levels)]

PREFIX = "# variance "
for line in stdin:
	if not line.startswith(PREFIX): continue
	level, var = line[len(PREFIX):].split()
	values[levels-int(level)].append(float(var))

level  = np.arange(levels) + logn
brown  = (np.power(2.0, 1-2*hurst) - 0.5) * np.power(2.0, -2*level*hurst);
length = np.fromiter(map(len, values), int, levels)
mean   = np.fromiter(map(np.mean, values), float, levels)
std    = np.fromiter(map(np.std, values), float, levels)

plt.figure(1)
plt.plot(level, brown, '-r')
plt.plot(level, mean, '+k')
plt.yscale('log')
plt.figure(2)
plt.errorbar(level, mean/brown, std/np.sqrt(length)/brown, fmt='xk')
plt.show()
