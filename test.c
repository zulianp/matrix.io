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
             "data/test/lhs.rowindex.raw",
             "data/test/lhs.colindex.raw",
             "data/test/lhs.value.raw",
             MPI_LONG,
             MPI_INT,
             MPI_FLOAT,
             &crs);

    crs_free(&crs);
    return MPI_Finalize();
}
