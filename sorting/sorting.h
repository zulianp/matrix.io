#ifndef MATRIXIO_SORTING_H
#define MATRIXIO_SORTING_H

#include <stddef.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

void local_argsort(ptrdiff_t size, MPI_Datatype key_type, void *const keys, ptrdiff_t *idx);
void gather_remap(ptrdiff_t size, int type_size, const ptrdiff_t *const gather_idx, void *const input, void *const output);

#ifdef __cplusplus
}
#endif

#endif  // MATRIXIO_SORTING_H
