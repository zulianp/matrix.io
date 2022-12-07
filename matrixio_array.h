#ifndef MATRIX_IO_ARRAY_H
#define MATRIX_IO_ARRAY_H

#include <mpi.h>
#include <stddef.h>

int array_read(MPI_Comm comm, const char *path, MPI_Datatype type, void **data, ptrdiff_t *out_nlocal, ptrdiff_t *out_nglobal);
int array_write(MPI_Comm comm, const char *path, MPI_Datatype type, const void *data, ptrdiff_t nlocal, ptrdiff_t out_nglobal);

#endif  // MATRIX_IO_ARRAY_H
