#!/usr/bin/env python3

import numpy as np; 
import sys

if len(sys.argv) != 6:
	print(f'usage {sys.argv[0]} <range_start> <range_end> <array_type> <array_file.array_type.raw> <output_file.array_type.raw>')
	exit(1)

range_start = np.int64(sys.argv[1])
range_end = np.int64(sys.argv[2])
array_type = sys.argv[3]
array_file = sys.argv[4]
output_file = sys.argv[5]

a = np.fromfile(array_file, dtype=np.dtype(array_type))
o = a[range_start:range_end]
o.tofile(output_file)
