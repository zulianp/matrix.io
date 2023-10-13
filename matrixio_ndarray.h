#ifndef MATRIX_IO_ND_ARRAY_H
#define MATRIX_IO_ND_ARRAY_H

#include <mpi.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int ndarray_read(MPI_Comm comm,
                 const char *path,
                 MPI_Datatype type,
                 int ndims,
                 void *data,
                 ptrdiff_t *const nlocal,
                 const ptrdiff_t *const nglobal);

int ndarray_write(MPI_Comm comm,
                  const char *path,
                  MPI_Datatype type,
                  int ndims,
                  const void *data,
                  const ptrdiff_t *const nlocal,
                  const ptrdiff_t *const nglobal);

int ndarray_read_segmented(MPI_Comm comm,
                           const char *path,
                           MPI_Datatype type,
                           int ndims,
                           void *data,
                           int segment_size,
                           ptrdiff_t *const nlocal,
                           const ptrdiff_t *const nglobal);

int ndarray_write_segmented(MPI_Comm comm,
                            const char *path,
                            MPI_Datatype type,
                            int ndims,
                            const void *data,
                            const int segment_size,
                            const ptrdiff_t *const nlocal,
                            const ptrdiff_t *const nglobal);

#ifdef __cplusplus
}
#endif

#endif  // MATRIX_IO_ND_ARRAY_H
