#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

    MPI_Barrier(comm);

    int rowptr_type_size;
    int colidx_type_size;
    int values_type_size;
    CATCH_MPI_ERROR(MPI_Type_size(crs.rowptr_type, &rowptr_type_size));
    CATCH_MPI_ERROR(MPI_Type_size(crs.colidx_type, &colidx_type_size));
    CATCH_MPI_ERROR(MPI_Type_size(crs.values_type, &values_type_size));

    ptrdiff_t nnz = crs.lnnz;

    for(int r = 0; r < size; r++) {

        if(r == rank) {
            printf("[%d] lnnz=%ld\n", rank, nnz);
        }

        fflush(stdout);
        MPI_Barrier(comm);
    }

    fflush(stdout);
    MPI_Barrier(comm);


    if (MATRIXIO_DENSE_OUTPUT) {
        
        ptrdiff_t maxcol = 0;
        for (ptrdiff_t k = 0; k < nnz; k++) {
            ptrdiff_t col = to_ptrdiff_t(colidx_type, &crs.colidx[k * colidx_type_size]);
            maxcol = MAX(col, maxcol);
        }

        maxcol += 1;
        double *buff = malloc(sizeof(double) * maxcol);

        for (int r = 0; r < size; ++r) {
            if (r == rank) {
                if (MATRIXIO_SHOW_PROCESS_RANK) {
                    printf("[%d]\n", rank);
                }

                for (ptrdiff_t i = 0; i < crs.lrows; ++i) {
                    memset(buff, 0, sizeof(double) * maxcol);
                    ptrdiff_t begin = to_ptrdiff_t(rowptr_type, &crs.rowptr[i * rowptr_type_size]) - crs.start;
                    ptrdiff_t end = to_ptrdiff_t(rowptr_type, &crs.rowptr[(i + 1) * rowptr_type_size]) - crs.start;

                    for (ptrdiff_t k = begin; k < end; k++) {
                        ptrdiff_t col = to_ptrdiff_t(colidx_type, &crs.colidx[k * colidx_type_size]);
                        double val = to_double(values_type, &crs.values[k * values_type_size]);
                        buff[col] = val;
                    }

                    for (ptrdiff_t c = 0; c < maxcol; c++) {
                        if(c + 1 < maxcol) {
                            printf("%g%s", buff[c], MATRIXIO_PRINT_CSV? "," : " ");
                        } else {
                            printf("%g", buff[c]);
                        }
                    }

                    printf("\n");
                }
            }

            fflush(stdout);
            MPI_Barrier(comm);
        }
    } else {
        for (int r = 0; r < size; ++r) {
            if (r == rank) {
                if (MATRIXIO_SHOW_PROCESS_RANK) {
                    printf("[%d]\n", rank);
                }
                for (ptrdiff_t i = 0; i < crs.lrows; ++i) {
                    ptrdiff_t begin = to_ptrdiff_t(rowptr_type, &crs.rowptr[i * rowptr_type_size]) - crs.start;
                    ptrdiff_t end = to_ptrdiff_t(rowptr_type, &crs.rowptr[(i + 1) * rowptr_type_size]) - crs.start;

                    printf("row %ld)\n", crs.rowoffset + i);
                    printf("cols: ");

                    for (ptrdiff_t k = begin; k < end; k++) {
                        printf("%ld ", to_ptrdiff_t(colidx_type, &crs.colidx[k * colidx_type_size]));
                    }

                    printf("\n");

                    printf("vals: ");
                    for (ptrdiff_t k = begin; k < end; k++) {
                        printf("%g ", to_double(values_type, &crs.values[k * values_type_size]));
                    }

                    printf("\n");
                }
                printf("\n");
            }

            MPI_Barrier(comm);
        }
    }

    crs_free(&crs);
    return MPI_Finalize();
}
