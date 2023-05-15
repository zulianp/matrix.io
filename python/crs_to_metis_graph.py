#!/usr/bin/env python3

import numpy as np
import scipy as sp
import sys
import os
import math

idx_t = np.int32
real_t = np.float64

def main(argv):
	if len(argv) != 3:
		print(f'usage: {argv[0]} <crs_folder> <out.graph>')
		exit(1)

	folder = argv[1]
	output_graph = argv[2]

	rowptr = np.fromfile(f'{folder}/rowptr.raw', dtype=idx_t)
	colidx = np.fromfile(f'{folder}/colidx.raw', dtype=idx_t)

	N = len(rowptr) - 1

	nedges = idx_t((len(colidx) - N) / 2)

	assert nedges * 2 == (len(colidx) - N)

	with open(output_graph, 'w') as f:	
		f.write(f'{N} {nedges}\n')
		for i in range(0, N):
			cols = colidx[rowptr[i]:rowptr[i+1]]

			for c in cols:
				if c != i:
					f.write(f'{c} ')

			f.write('\n')
		f.write('# End file')

if __name__ == '__main__':
	main(sys.argv)
