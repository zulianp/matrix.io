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
    ptrdiff_t lnnz;
    ptrdiff_t gnnz;
    ptrdiff_t start;
    ptrdiff_t rowoffset;

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

int crs_read_str(MPI_Comm comm,
                 const char *rowptr_path,
                 const char *colidx_path,
                 const char *values_path,
                 const char *rowptr_type,
                 const char *colidx_type,
                 const char *values_type,
                 crs_t *crs);

int crs_write(MPI_Comm comm,
              const char *rowptr_path,
              const char *colidx_path,
              const char *values_path,
              MPI_Datatype rowptr_type,
              MPI_Datatype colidx_type,
              MPI_Datatype values_type,
              crs_t *crs);

int crs_read_folder(MPI_Comm comm,
              const char *folder,
              MPI_Datatype rowptr_type,
              MPI_Datatype colidx_type,
              MPI_Datatype values_type,
              crs_t *crs);

int crs_write_folder(MPI_Comm comm,
              const char *folder,
              MPI_Datatype rowptr_type,
              MPI_Datatype colidx_type,
              MPI_Datatype values_type,
              crs_t *crs);

int crs_write_str(MPI_Comm comm,
                  const char *rowptr_path,
                  const char *colidx_path,
                  const char *values_path,
                  const char *rowptr_type,
                  const char *colidx_type,
                  const char *values_type,
                  crs_t *crs);

/// Free memory
int crs_free(crs_t *crs);

// Memory is managed outside
int crs_release(crs_t *crs);

#endif  // MATRIX_IO_CRS_H
