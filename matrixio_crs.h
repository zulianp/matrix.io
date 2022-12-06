#ifndef MATRIX_IO_CRS_H
#define MATRIX_IO_CRS_H

#include <mpi.h>
#include <stddef.h>

typedef struct {
    char *rowptr;
    char *colidx;
    char *values;

    ptrdiff_t grows;
    ptrdiff_t lrows;
    ptrdiff_t nnz;
    ptrdiff_t start;

    int rowptr_type_size;
    int colidx_type_size;
    int values_type_size;
} crs_t;

int crs_read(MPI_Comm comm,
             const char *rowptr_path,
             const char *colidx_path,
             const char *values_path,
             MPI_Datatype rowptr_type,
             MPI_Datatype colidx_type,
             MPI_Datatype values_type,
             crs_t *crs);

/// Free memory
int crs_free(crs_t *crs);

// Memory is managed outside
int crs_release(crs_t *crs);

#endif  // MATRIX_IO_CRS_H
