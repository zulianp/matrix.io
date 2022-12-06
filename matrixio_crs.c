#include "matrixio_crs.h"

#include "utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int crs_read_str(MPI_Comm comm,
                 const char *rowptr_path,
                 const char *colidx_path,
                 const char *values_path,
                 const char *rowptr_type_str,
                 const char *colidx_type_str,
                 const char *values_type_str,
                 crs_t *crs) {
    MPI_Datatype rowptr_type = string_to_mpi_datatype(rowptr_type_str);
    MPI_Datatype colidx_type = string_to_mpi_datatype(colidx_type_str);
    MPI_Datatype values_type = string_to_mpi_datatype(values_type_str);
    return crs_read(comm, rowptr_path, colidx_path, values_path, rowptr_type, colidx_type, values_type, crs);
}

int crs_read(MPI_Comm comm,
             const char *rowptr_path,
             const char *colidx_path,
             const char *values_path,
             MPI_Datatype rowptr_type,
             MPI_Datatype colidx_type,
             MPI_Datatype values_type,
             crs_t *crs) {
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_File file;
    MPI_Offset rowptr_nbytes = -1;
    int rowptr_type_size = 0;

    CATCH_MPI_ERROR(MPI_File_open(comm, rowptr_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    CATCH_MPI_ERROR(MPI_File_get_size(file, &rowptr_nbytes));
    CATCH_MPI_ERROR(MPI_Type_size(rowptr_type, &rowptr_type_size));

    ptrdiff_t nrows = rowptr_nbytes / rowptr_type_size;

    if (nrows * rowptr_type_size != rowptr_nbytes) {
        if (!rank) {
            fprintf(stderr, "[Error] Wrong type specification for rowptr\n");
        }

        MPI_Abort(comm, 1);
    }

    // Remove last rowptr
    nrows -= 1;
    crs->grows = nrows;

    ptrdiff_t uniform_split = nrows / size;
    ptrdiff_t nlocal = uniform_split;
    ptrdiff_t remainder = nrows - nlocal * size;

    if (remainder > rank) {
        nlocal += 1;
    }

    ///////////////////////////////////////////////////////

    ptrdiff_t offset = rank * uniform_split;
    offset += MIN(rank, remainder);

    char *rowptr = (char *)malloc((nlocal + 1) * rowptr_type_size);

    MPI_Status status;

    ///////////////////////////////////////////////////////
    // Read rowptr
    ///////////////////////////////////////////////////////

    CATCH_MPI_ERROR(MPI_File_read_at_all(file, offset * rowptr_type_size, rowptr, nlocal + 1, rowptr_type, &status));
    CATCH_MPI_ERROR(MPI_File_close(&file));

    ///////////////////////////////////////////////////////
    // Read colidx
    ///////////////////////////////////////////////////////

    ptrdiff_t start = to_ptrdiff_t(rowptr_type, &rowptr[0]);
    ptrdiff_t end = to_ptrdiff_t(rowptr_type, &(rowptr[nlocal * rowptr_type_size]));

    crs->start = start;

    ptrdiff_t nnz = end - start;
    int colidx_type_size = 0;

    crs->nnz = nnz;
    CATCH_MPI_ERROR(MPI_Type_size(colidx_type, &colidx_type_size));

    char *colidx = (char *)malloc(nnz * colidx_type_size);

    CATCH_MPI_ERROR(MPI_File_open(comm, colidx_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));

    CATCH_MPI_ERROR(MPI_File_read_at_all(file, start * colidx_type_size, colidx, nnz, colidx_type, &status));

    CATCH_MPI_ERROR(MPI_File_close(&file));

    ///////////////////////////////////////////////////////
    // Read values
    ///////////////////////////////////////////////////////

    int values_type_size = 0;
    CATCH_MPI_ERROR(MPI_Type_size(values_type, &values_type_size));
    char *values = (char *)malloc(nnz * values_type_size);

    CATCH_MPI_ERROR(MPI_File_open(comm, values_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));

    CATCH_MPI_ERROR(MPI_File_read_at_all(file, start * values_type_size, values, nnz, values_type, &status));

    CATCH_MPI_ERROR(MPI_File_close(&file));

    ///////////////////////////////////////////////////////

    MPI_Barrier(comm);

    crs->rowptr = rowptr;
    crs->colidx = colidx;
    crs->values = values;

    crs->grows = nrows;
    crs->lrows = nlocal;
    crs->nnz = nnz;
    crs->start = start;

    crs->rowptr_type_size = rowptr_type_size;
    crs->colidx_type_size = colidx_type_size;
    crs->values_type_size = values_type_size;
    return 0;
}

int crs_free(crs_t *crs) {
    free(crs->rowptr);
    free(crs->colidx);
    free(crs->values);
    crs->lrows = 0;
    crs->grows = 0;
    crs->nnz = 0;
    crs->start = 0;
    return 0;
}

int crs_release(crs_t *crs) {
    crs->rowptr = 0;
    crs->colidx = 0;
    crs->values = 0;
    crs->lrows = 0;
    crs->grows = 0;
    crs->nnz = 0;
    crs->start = 0;
    return 0;
}
