#ifndef MATRIXIO_PARMETIS_H
#define MATRIXIO_PARMETIS_H

#include "matrixio_crs.h"

#include <mpi.h>

int crs_decompose(MPI_Comm comm, crs_t*crs, const int n_parts, int*const parts);
int crs_graph_decompose(MPI_Comm comm, crs_graph_t*crs, const int n_parts, int*const parts);

#endif //MATRIXIO_PARMETIS_H
