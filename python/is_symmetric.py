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
		# print(argv)
		print(f'usage: {argv[0]} <crs_folder> <check_data>')
		exit(1)

	folder = argv[1]
	check_data = int(argv[2])

	rowptr = np.fromfile(f'{folder}/rowptr.raw', dtype=idx_t)
	colidx = np.fromfile(f'{folder}/colidx.raw', dtype=idx_t)

	N = len(rowptr) - 1

	if check_data:
		data   = np.fromfile(f'{folder}/values.raw', dtype=real_t)
	else:
		data = np.ones(colidx.shape, dtype=real_t)

	A = sp.sparse.csr_matrix((data, colidx, rowptr), shape=(N, N)) 
	D = A.T - A

	nnz = np.sum(np.abs(D.data))

	print(f'sum(abs(A.T - A)) = {nnz}')
	# print(D.data)

if __name__ == '__main__':
	main(sys.argv)
