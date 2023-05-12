#ifndef MATRIXIO_PARMETIS_H
#define MATRIXIO_PARMETIS_H

#include "matrixio_crs.h"

#include <mpi.h>

int decompose(MPI_Comm comm, crs_t*crs, const int n_parts, int*const parts);

#endif //MATRIXIO_PARMETIS_H
