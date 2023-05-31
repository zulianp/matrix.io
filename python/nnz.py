#!/usr/bin/env python3

import numpy as np
import scipy as sp
import sys
import os
import math

idx_t = np.int32
real_t = np.float64

def main(argv):
	if len(argv) != 2:
		print(f'usage: {argv[0]} <data.raw>')
		exit(1)

	path = argv[1]
	data = np.fromfile(path, dtype=real_t)
	print(f'nnz={np.count_nonzero(data)}/{len(data)}')

if '__main__' == __name__:
	main(sys.argv)