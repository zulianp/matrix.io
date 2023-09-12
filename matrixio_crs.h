#ifndef MATRIX_IO_CRS_H
#define MATRIX_IO_CRS_H

#include <mpi.h>
#include <stddef.h>

#include "matrixio_base.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    matrixio_byte_t *rowptr;
    matrixio_byte_t *colidx;

    ptrdiff_t grows;
    ptrdiff_t lrows;
    ptrdiff_t lnnz;
    ptrdiff_t gnnz;
    ptrdiff_t start;
    ptrdiff_t rowoffset;

    MPI_Datatype rowptr_type;
    MPI_Datatype colidx_type;
} crs_graph_t;

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

    MPI_Datatype rowptr_type;
    MPI_Datatype colidx_type;
    MPI_Datatype values_type;
} crs_t;

typedef struct {
    char *rowptr;
    char *colidx;

    int block_size;
    // SoA values
    char **values;

    ptrdiff_t grows;
    ptrdiff_t lrows;
    ptrdiff_t lnnz;
    ptrdiff_t gnnz;
    ptrdiff_t start;
    ptrdiff_t rowoffset;

    MPI_Datatype rowptr_type;
    MPI_Datatype colidx_type;
    MPI_Datatype values_type;
} block_crs_t;

int crs_alloc_same(const crs_t *const tpl, crs_t *result);

// graph is rendered invalid due to pointer ownership passed to crs
int crs_graph_to_crs(crs_graph_t *const graph, crs_t *const crs);
int crs_graph_to_block_crs(crs_graph_t *const graph, block_crs_t *const crs);

int crs_graph_view_from_crs(crs_t *const crs, crs_graph_t *const graph);
int crs_graph_view_from_block_crs(block_crs_t *const crs, crs_graph_t *const graph);

int crs_graph_read(MPI_Comm comm,
                   const char *rowptr_path,
                   const char *colidx_path,
                   MPI_Datatype rowptr_type,
                   MPI_Datatype colidx_type,
                   crs_graph_t *crs);

int crs_graph_read_values(MPI_Comm comm,
                          const crs_graph_t *const crs,
                          const char *values_path,
                          MPI_Datatype values_type,
                          matrixio_byte_t *const values);

int crs_read(MPI_Comm comm,
             const char *rowptr_path,
             const char *colidx_path,
             const char *values_path,
             MPI_Datatype rowptr_type,
             MPI_Datatype colidx_type,
             MPI_Datatype values_type,
             crs_t *crs);

int crs_read_AoS_block(MPI_Comm comm,
                       const char *rowptr_path,
                       const char *colidx_path,
                       const char *values_path,
                       MPI_Datatype rowptr_type,
                       MPI_Datatype colidx_type,
                       MPI_Datatype values_type,
                       const int block_size,
                       crs_t *crs);

int block_crs_read(MPI_Comm comm,
                   const char *rowptr_path,
                   const char *colidx_path,
                   const char *values_pattern,
                   MPI_Datatype rowptr_type,
                   MPI_Datatype colidx_type,
                   MPI_Datatype values_type,
                   block_crs_t *crs);



int crs_read_str(MPI_Comm comm,
                 const char *rowptr_path,
                 const char *colidx_path,
                 const char *values_path,
                 const char *rowptr_type,
                 const char *colidx_type,
                 const char *values_type,
                 crs_t *crs);

int crs_read_folder(MPI_Comm comm,
                    const char *folder,
                    MPI_Datatype rowptr_type,
                    MPI_Datatype colidx_type,
                    MPI_Datatype values_type,
                    crs_t *crs);

int crs_write(MPI_Comm comm, const char *rowptr_path, const char *colidx_path, const char *values_path, crs_t *crs);
int block_crs_write(MPI_Comm comm, const char *rowptr_path, const char *colidx_path, const char *values_format, block_crs_t *crs);

int crs_write_folder(MPI_Comm comm, const char *folder, crs_t *crs);

int crs_graph_write(MPI_Comm comm, const char *rowptr_path, const char *colidx_path, crs_graph_t *crs);
int crs_graph_write_values(MPI_Comm comm,
                           const crs_graph_t *const crs,
                           const char *values_path,
                           MPI_Datatype values_type,
                           matrixio_byte_t *const values);

/// Free memory
int crs_free(crs_t *const crs);
int block_crs_free(block_crs_t *const crs);

// Memory is managed outside
int crs_release(crs_t *const crs);

int crs_graph_free(crs_graph_t *const crs);
int crs_graph_release(crs_graph_t *const crs);

#ifdef __cplusplus
}
#endif

#endif  // MATRIX_IO_CRS_H
