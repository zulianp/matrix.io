#ifndef MATRIX_IO_ARRAY_H
#define MATRIX_IO_ARRAY_H

#include <mpi.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int array_create_from_file_segmented(MPI_Comm comm,
                                     const char *path,
                                     MPI_Datatype type,
                                     void **data,
                                     const int segment_size,
                                     ptrdiff_t *out_nlocal,
                                     ptrdiff_t *out_nglobal);

int array_create_from_file(MPI_Comm comm,
                           const char *path,
                           MPI_Datatype type,
                           void **data,
                           ptrdiff_t *out_nlocal,
                           ptrdiff_t *out_nglobal);

int array_read(MPI_Comm comm, const char *path, MPI_Datatype type, void *data, ptrdiff_t nlocal, ptrdiff_t nglobal);

int array_read_segmented(MPI_Comm comm,
                         const char *path,
                         MPI_Datatype type,
                         void *data,
                         int segment_size,
                         ptrdiff_t nlocal,
                         ptrdiff_t nglobal);

int array_write(MPI_Comm comm,
                const char *path,
                MPI_Datatype type,
                const void *data,
                ptrdiff_t nlocal,
                ptrdiff_t out_nglobal);

int array_write_segmented(MPI_Comm comm,
                          const char *path,
                          MPI_Datatype type,
                          const void *data,
                          const int segment_size,
                          ptrdiff_t nlocal,
                          ptrdiff_t nglobal);

int array_range_select(MPI_Comm comm,
                       MPI_Datatype type,
                       void *const in,
                       void *const out,
                       ptrdiff_t in_nlocal,
                       ptrdiff_t range_start,
                       ptrdiff_t range_end);

int array_convert(
    const ptrdiff_t size, 
    MPI_Datatype from_type,
    const void * from, 
    MPI_Datatype to_type, 
    void * to);

#ifdef __cplusplus
}
#endif

#endif  // MATRIX_IO_ARRAY_H
