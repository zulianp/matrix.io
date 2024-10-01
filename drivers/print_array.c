#include <stdio.h>
#include <stdlib.h>

#include "matrixio_array.h"
#include "matrixio_crs.h"

#include "utils.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (argc < 2) {
        if (!rank) {
            fprintf(stderr, "usage: %s <array.raw> [type=float]", argv[0]);
        }

        return EXIT_FAILURE;
    }

    MPI_Datatype type = MPI_FLOAT;

    if (argc > 2) {
        type = string_to_mpi_datatype(argv[2]);
    }


    int MATRIXIO_SHOW_PROCESS_RANK = 1;
    MATRIXIO_READ_ENV(MATRIXIO_SHOW_PROCESS_RANK, atoi);

    int MATRIXIO_PRINT_CSV = 0;
    MATRIXIO_READ_ENV(MATRIXIO_PRINT_CSV, atoi);

    ptrdiff_t max_entries_x_line = 50;
    ptrdiff_t nlocal, nglobal;
    char *data;
    array_create_from_file(comm, argv[1], type, (void **)&data, &nlocal, &nglobal);

    int type_size;
    MPI_CATCH_ERROR(MPI_Type_size(type, &type_size));

    MPI_Barrier(comm);

    for (int r = 0; r < size; ++r) {
        if (r == rank) {
            if(MATRIXIO_SHOW_PROCESS_RANK) {
                printf("[%d]\n", rank);
            }

            for (ptrdiff_t i = 0; i < nlocal; ++i) {
                double v = to_double(type, &data[i * type_size]);

                if(MATRIXIO_PRINT_CSV) {
                    if(i == nlocal - 1) {
                        printf("%g ", v);
                    }else {
                        printf("%g, ", v);
                    }
                } else {
                    if((i+1) % max_entries_x_line == 0) {
                        printf("\n");
                    }

                    printf("%g ", v);
                }


            }
            printf("\n");
        }

        MPI_Barrier(comm);
    }

    free(data);
    return MPI_Finalize();
}
