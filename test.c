#include <stdio.h>
#include "matrixio_crs.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    crs_t crs;
    crs_read(comm,
             "/Users/patrickzulian/Desktop/code/splinsol/Sources/solve/test/"
             "lhs.rowindex.raw",
             "/Users/patrickzulian/Desktop/code/splinsol/Sources/solve/test/"
             "lhs.colindex.raw",
             "/Users/patrickzulian/Desktop/code/splinsol/Sources/solve/test/"
             "lhs.value.raw",
             MPI_LONG,
             MPI_INT,
             MPI_FLOAT,
             &crs);

    long *rowptr = (long *)crs.rowptr;
    int *colidx = (int *)crs.colidx;

    MPI_Barrier(comm);

    for (int r = 0; r < size; ++r) {
        if (r == rank) {
            printf("[%d] ", rank);
            for (ptrdiff_t i = 0; i < crs.lrows; ++i) {
                printf("%ld) %ld-%ld: ", i, rowptr[i], rowptr[i + 1]);

                for (long k = rowptr[i]; k < rowptr[i + 1]; k++) {
                    long idx = k - crs.start;
                    printf("%d ", colidx[idx]);
                }
                printf("\n");
            }
            printf("\n");
        }

        MPI_Barrier(comm);
    }

    crs_free(&crs);
    return MPI_Finalize();
}
