#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "matrixio_array.h"
#include "matrixio_crs.h"
#include "matrixio_parmetis.h"

#include "matrixio_checks.h"

#include "utils.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (argc < 3 || argc > 5) {
        if (!rank) {
            fprintf(stderr,
                    "usage: %s <rowptr.raw> <colidx.raw> [rowptr_type=int] [colidx_type=int]",
                    argv[0]);
        }

        return EXIT_FAILURE;
    }

    int MATRIXIO_NPARTS = size;
    MATRIXIO_READ_ENV(MATRIXIO_NPARTS, atoi);

    MPI_Datatype rowptr_type = MPI_INT;
    MPI_Datatype colidx_type = MPI_INT;

    if (argc > 3) {
        rowptr_type = string_to_mpi_datatype(argv[3]);
    }

    if (argc > 4) {
        colidx_type = string_to_mpi_datatype(argv[4]);
    }

    crs_graph_t crs;
    crs_graph_read(comm, argv[1], argv[2], rowptr_type, colidx_type, &crs);

    int err = crs_check_graph(crs.lrows, crs.grows, crs.lnnz, crs.rowptr_type, crs.rowptr, crs.colidx_type, crs.colidx);
    MPI_Allreduce(MPI_IN_PLACE, &err, 1, MPI_INT, MPI_MAX, comm);

    if (err) {
        if (!rank) {
            printf("crs_check_graph returned error code %d. Aborting...\n", err);
        }

        MPI_Abort(comm, 1);
    }

    int *parts = (int *)malloc(crs.lrows * sizeof(int));
    crs_graph_decompose(comm, &crs, MATRIXIO_NPARTS, parts);
    array_write(comm, "parts.int32.raw", MPI_INT, parts, crs.lrows, crs.grows);

    crs_graph_free(&crs);
    return MPI_Finalize();
}
