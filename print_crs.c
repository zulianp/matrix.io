#include <stdio.h>
#include <stdlib.h>

#include "matrixio_crs.h"

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

    MPI_Barrier(comm);

    for (int r = 0; r < size; ++r) {
        if (r == rank) {
            printf("[%d]\n", rank);
            for (ptrdiff_t i = 0; i < crs.lrows; ++i) {
                ptrdiff_t begin = to_ptrdiff_t(rowptr_type, &crs.rowptr[i * crs.rowptr_type_size]) - crs.start;
                ptrdiff_t end = to_ptrdiff_t(rowptr_type, &crs.rowptr[(i + 1) * crs.rowptr_type_size]) - crs.start;

                printf("row %ld)\n", i);
                printf("cols: ");

                for (ptrdiff_t k = begin; k < end; k++) {
                    printf("%ld ", to_ptrdiff_t(colidx_type, &crs.colidx[k * crs.colidx_type_size]));
                }

                printf("\n");

                printf("vals: ");
                for (ptrdiff_t k = begin; k < end; k++) {
                    printf("%g ", to_double(values_type, &crs.values[k * crs.values_type_size]));
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
