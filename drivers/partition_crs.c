#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrixio_crs.h"
#include "matrixio_array.h"
#include "matrixio_parmetis.h"

#include "utils.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (argc < 4) {
        if (!rank) {
            fprintf(stderr,
                    "usage: %s <rowptr.raw> <colidx.raw> <values.raw> [rowptr_type=int] [colidx_type=int] "
                    "[values_type=float]",
                    argv[0]);
        }

        return EXIT_FAILURE;
    }

    int MATRIXIO_DENSE_OUTPUT = 0;
    MATRIXIO_READ_ENV(MATRIXIO_DENSE_OUTPUT, atoi);

    int MATRIXIO_SHOW_PROCESS_RANK = 1;
    MATRIXIO_READ_ENV(MATRIXIO_SHOW_PROCESS_RANK, atoi);

    int MATRIXIO_PRINT_CSV = 0;
    MATRIXIO_READ_ENV(MATRIXIO_PRINT_CSV, atoi);

    MPI_Datatype rowptr_type = MPI_INT;
    MPI_Datatype colidx_type = MPI_INT;
    MPI_Datatype values_type = MPI_FLOAT;

    if (argc > 4) {
        rowptr_type = string_to_mpi_datatype(argv[4]);
    }

    if (argc > 5) {
        colidx_type = string_to_mpi_datatype(argv[5]);
    }

    if (argc > 6) {
        values_type = string_to_mpi_datatype(argv[6]);
    }

    crs_t crs;
    crs_read(comm, argv[1], argv[2], argv[3], rowptr_type, colidx_type, values_type, &crs);

    int *parts = (int *)malloc(crs.lrows * sizeof(int));
    decompose(comm, &crs, parts);

    array_write(comm, "parts.int32.raw", MPI_INT, parts, crs.lrows, crs.grows);

    crs_free(&crs);
    return MPI_Finalize();
}
