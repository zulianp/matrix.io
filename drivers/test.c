#include "matrixio_array.h"
#include "matrixio_crs.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    {
        // Run read crs
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
                  &crs);

        crs_free(&crs);
    }

    {
        // Run read vector
        ptrdiff_t nlocal, nglobal;
        char *data;
        array_create_from_file(comm, "data/test/rhs.raw", MPI_FLOAT, (void **)&data, &nlocal, &nglobal);
        array_write(comm, "data/test/dump.raw", MPI_FLOAT, data, nlocal, nglobal);
        array_write_segmented(comm, "data/test/dump_seg.raw", MPI_FLOAT, data, 2, nlocal, nglobal);
        free(data);
    }

    {
        ptrdiff_t nlocal = 10;
        ptrdiff_t nglobal = nlocal * size;
        double *data = malloc(nlocal * sizeof(double));
        float *test_data = malloc(nlocal * sizeof(float));

        for(ptrdiff_t i = 0; i < nlocal; i++) {
            test_data[i] = i;
        }

        array_write_convert(comm, "data/test/dump.float32", MPI_FLOAT, test_data, nlocal, nglobal);

        // Read and convert to MPI type (float32 -> double)
        array_read_convert(comm,  "data/test/dump.float32", MPI_DOUBLE, data, nlocal, nglobal);
        array_write_convert(comm, "data/test/dump.float32", MPI_DOUBLE, data, nlocal, nglobal);
        array_read_convert(comm,  "data/test/dump.float32", MPI_DOUBLE, data, nlocal, nglobal);

        for(ptrdiff_t i = 0; i < nlocal; i++) {
            int expected = test_data[i];
            int actual = data[i];

            assert(expected == actual);
            if(expected != actual) {
                printf("%f != %g\n", test_data[i], data[i]);
            }
        }

        free(data);
        free(test_data);
    }

    return MPI_Finalize();
}
