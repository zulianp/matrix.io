#include "matrixio_array.h"
#include "matrixio_crs.h"

#include <stdio.h>
#include <stdlib.h> 

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    {
        //Run read crs
        crs_t crs;
        crs_read(comm,
                 "data/test/lhs.rowindex.raw",
                 "data/test/lhs.colindex.raw",
                 "data/test/lhs.value.raw",
                 MPI_LONG,
                 MPI_INT,
                 MPI_FLOAT,
                 &crs);


        crs_write(comm,
                 "data/test/dump.rowindex.raw",
                 "data/test/dump.colindex.raw",
                 "data/test/dump.value.raw",
                 MPI_LONG,
                 MPI_INT,
                 MPI_FLOAT,
                 &crs);

        crs_free(&crs);
    }

    {
        //Run read rhs
        ptrdiff_t nlocal, nglobal;
        char *data;
        array_read(comm, "data/test/rhs.raw", MPI_FLOAT, (void**)&data, &nlocal, &nglobal);

        array_write(comm, "data/test/dump.raw", MPI_FLOAT, data, nlocal, nglobal);

        free(data);
    }

    return MPI_Finalize();
}
