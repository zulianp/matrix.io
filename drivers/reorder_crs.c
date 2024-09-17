#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "matrixio_array.h"
#include "matrixio_crs.h"
#include "matrixio_parmetis.h"

#include "matrixio_checks.h"

#include "sorting.h"
#include "utils.h"

static int array_max_int(const int *array, ptrdiff_t n) {
    if (!n) return 0;

    int ret = array[0];
    for (ptrdiff_t i = 1; i < n; i++) {
        ret = MAX(ret, array[i]);
    }

    return ret;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (argc < 4) {
        if (!rank) {
            fprintf(stderr,
                    "usage: %s <parts.raw> <rowptr.raw> <colidx.raw> <values.raw> [rowptr_type=int] "
                    "[colidx_type=int] "
                    "[values_type=float]",
                    argv[0]);
        }

        return EXIT_FAILURE;
    }

    if (size != 1) {
        fprintf(stderr, "Only serial runs are supported for now!\n");
        return EXIT_FAILURE;
    }

    int MATRIXIO_NPARTS = size;
    MATRIXIO_READ_ENV(MATRIXIO_NPARTS, atoi);

    int MATRIXIO_SORT_COLUMNS = 1;
    MATRIXIO_READ_ENV(MATRIXIO_SORT_COLUMNS, atoi);

    // FIXME
    typedef int rowptr_t;
    typedef int colidx_t;

    MPI_Datatype parts_type = MPI_INT;
    MPI_Datatype rowptr_type = MPI_INT;
    MPI_Datatype colidx_type = MPI_INT;
    MPI_Datatype values_type = MPI_FLOAT;

    if (argc > 5) {
        rowptr_type = string_to_mpi_datatype(argv[5]);
    }

    if (argc > 6) {
        colidx_type = string_to_mpi_datatype(argv[6]);
    }

    if (argc > 7) {
        values_type = string_to_mpi_datatype(argv[7]);
    }

    int values_type_size = 0;
    MPI_CATCH_ERROR(MPI_Type_size(values_type, &values_type_size));

    int colidx_type_size = 0;
    MPI_CATCH_ERROR(MPI_Type_size(colidx_type, &colidx_type_size));

    ptrdiff_t nl = 0, ng = 0;
    int *parts = 0;
    array_create_from_file(comm, argv[1], parts_type, (void **)&parts, &nl, &ng);

    int max_part_id = array_max_int(parts, nl);
    ptrdiff_t nn = MAX(max_part_id, size);

    ptrdiff_t *part_ptr = (ptrdiff_t *)calloc(nn + 1, sizeof(ptrdiff_t));
    ptrdiff_t *count = (ptrdiff_t *)calloc(nn + 1, sizeof(ptrdiff_t));

    for (ptrdiff_t i = 0; i < nl; i++) {
        part_ptr[parts[i] + 1]++;
    }

    for (ptrdiff_t i = 0; i < nn; i++) {
        part_ptr[i + 1] += part_ptr[i];
    }

    crs_t crs;
    crs_read(comm, argv[2], argv[3], argv[4], rowptr_type, colidx_type, values_type, &crs);

    ptrdiff_t *mapping = (ptrdiff_t *)malloc(nl * sizeof(ptrdiff_t));

    // array_write(comm, "parts.int32.raw", MPI_INT, parts, crs.lrows, crs.grows);

    // Create order preserving mapping
    for (ptrdiff_t i = 0; i < nl; i++) {
        const int p = parts[i];
        ptrdiff_t begin = part_ptr[p];
        mapping[i] = begin + count[p]++;
    }

    crs_t result;
    crs_alloc_same(&crs, &result);

    const rowptr_t *const crs_rowptr = (rowptr_t *)crs.rowptr;
    rowptr_t *const result_rowptr = (rowptr_t *)result.rowptr;

    const colidx_t *const crs_colidx = (colidx_t *)crs.colidx;
    colidx_t *const result_colidx = (colidx_t *)result.colidx;

    result_rowptr[0] = 0;
    for (ptrdiff_t i = 0; i < nl; i++) {
        result_rowptr[mapping[i] + 1] = crs_rowptr[i + 1] - crs_rowptr[i];
    }

    for (ptrdiff_t i = 0; i < nl; i++) {
        result_rowptr[i + 1] += result_rowptr[i];
    }

    rowptr_t max_row_len = 0;
    for (ptrdiff_t i = 0; i < nl; i++) {
        const rowptr_t crs_begin = crs_rowptr[i];
        const rowptr_t crs_extent = crs_rowptr[i + 1] - crs_begin;

        const rowptr_t result_begin = result_rowptr[i];
        const rowptr_t result_extent = result_rowptr[i + 1] - result_begin;

        const colidx_t *crs_row = &crs_colidx[crs_begin];
        colidx_t *result_row = &result_colidx[result_begin];

        assert(crs_extent == result_extent);

        max_row_len = MAX(max_row_len, crs_extent);

        for (rowptr_t i = 0; i < crs_extent; i++) {
            result_row[i] = mapping[crs_row[i]];
        }

        for (rowptr_t i = 0; i < crs_extent; i++) {
            memcpy(&result.values[i * values_type_size], &crs.values[i * values_type_size], values_type_size);
        }
    }

    if (MATRIXIO_SORT_COLUMNS) {
        ptrdiff_t *idx = malloc(max_row_len * sizeof(ptrdiff_t));
        uint8_t *buff = malloc(max_row_len * MAX(values_type_size, colidx_type_size));

        for (ptrdiff_t i = 0; i < nl; i++) {
            const rowptr_t result_begin = result_rowptr[i];
            const rowptr_t result_extent = result_rowptr[i + 1] - result_begin;
            colidx_t *result_row = &result_colidx[result_begin];

            // Compute reordering index
            local_argsort(result_extent, colidx_type, result_row, idx);

            memcpy(buff, result_row, result_extent * colidx_type_size);
            gather_remap(result_extent, colidx_type_size, idx, buff, result_row);

            memcpy(buff, &result.values[result_begin * values_type_size], result_extent * values_type_size);
            gather_remap(result_extent, values_type_size, idx, buff, &result.values[result_begin * values_type_size]);
        }

        free(idx);
        free(buff);
    }

    array_write(comm, "reorder_mapping.raw", MPI_LONG, mapping, crs.lrows, crs.grows);

    // FIXME
    array_write(comm, "parallel_layout.raw", MPI_LONG, part_ptr, size+1, size+1);

    {
        // Clean-up
        crs_free(&crs);
        crs_free(&result);
        free(parts);
        free(part_ptr);
        free(count);
        free(mapping);
    }

    return MPI_Finalize();
}
