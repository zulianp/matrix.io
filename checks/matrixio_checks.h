#ifndef MATRIXIO_CHECK_H
#define MATRIXIO_CHECK_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <mpi.h>


/**
 * @brief Checks the integrity of a compressed row storage (CRS) graph represented by its row pointer and column index arrays.
 *
 * This function verifies the correctness of the row pointer array and column index array in terms of their values and sizes.
 *
 * @param nrows The number of rows in the CRS graph.
 * @param ncols The number of columns in the CRS graph.
 * @param nnz The number of non-zero elements in the CRS graph.
 * @param rowptr_type MPI datatype of the row pointer array (MPI_INT, MPI_LONG, or MPI_UNSIGNED).
 * @param rowptr Pointer to the row pointer array.
 * @param colidx_type MPI datatype of the column index array (MPI_INT, MPI_LONG, or MPI_UNSIGNED).
 * @param colidx Pointer to the column index array.
 *
 * @return An integer code indicating the result of the check:
 *     - 0: The CRS graph is sane, i.e., the row pointer and column index arrays are correct.
 *     - 1: The row pointer array is not sane.
 *     - 2: The column index array is not sane.
 * 	   - 3: Both are not sane.  	
 */  
int crs_check_graph(const ptrdiff_t nrows,
                    const ptrdiff_t ncols,
                    const ptrdiff_t nnz,
                    MPI_Datatype rowptr_type,
                    const void *const rowptr,
                    MPI_Datatype colidx_type,
                    const void *const colidx);


#ifdef __cplusplus
}
#endif

#endif //MATRIXIO_CHECK_H

