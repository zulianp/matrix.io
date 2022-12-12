#include "matrixio_crs.h"

#include "utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

    MPI_Offset gnnz_bytes = -1;
    ptrdiff_t start = to_ptrdiff_t(rowptr_type, &rowptr[0]);
    ptrdiff_t end = to_ptrdiff_t(rowptr_type, &(rowptr[nlocal * rowptr_type_size]));

    crs->start = start;

    ptrdiff_t nnz = end - start;
    int colidx_type_size = 0;

    crs->lnnz = nnz;
    CATCH_MPI_ERROR(MPI_Type_size(colidx_type, &colidx_type_size));

    char *colidx = (char *)malloc(nnz * colidx_type_size);

    CATCH_MPI_ERROR(MPI_File_open(comm, colidx_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));

    CATCH_MPI_ERROR(MPI_File_get_size(file, &gnnz_bytes));

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

    // MPI_Barrier(comm);

    crs->rowptr = rowptr;
    crs->colidx = colidx;
    crs->values = values;

    crs->grows = nrows;
    crs->lrows = nlocal;
    crs->gnnz = gnnz_bytes / colidx_type_size;
    crs->lnnz = nnz;
    crs->start = start;

    // crs->rowptr_type_size = rowptr_type_size;
    // crs->colidx_type_size = colidx_type_size;
    // crs->values_type_size = values_type_size;

    crs->rowptr_type = rowptr_type;
    crs->colidx_type = colidx_type;
    crs->values_type = values_type;
    crs->rowoffset = offset;

    // printf("[read] grows=%ld nrows=%ld nlocal=%ld\n", (long)crs->grows, (long)nrows, (long)nlocal);
    return 0;
}

int crs_read_folder(MPI_Comm comm,
                    const char *folder,
                    MPI_Datatype rowptr_type,
                    MPI_Datatype colidx_type,
                    MPI_Datatype values_type,
                    crs_t *crs) {
    assert(strlen(folder) + 11 < 1024);

    char rowptr_path[1024];
    char colidx_path[1024];
    char values_path[1024];

    sprintf(rowptr_path, "%s/rowptr.raw", folder);
    sprintf(colidx_path, "%s/colidx.raw", folder);
    sprintf(values_path, "%s/values.raw", folder);

    return crs_read(comm, rowptr_path, colidx_path, values_path, rowptr_type, colidx_type, values_type, crs);
}

int crs_write(MPI_Comm comm, const char *rowptr_path, const char *colidx_path, const char *values_path, crs_t *crs) {
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_File file;
    MPI_Status status;

    MPI_Datatype rowptr_type = crs->rowptr_type;
    MPI_Datatype colidx_type = crs->colidx_type;
    MPI_Datatype values_type = crs->values_type;

    int rowptr_type_size = 0;
    int colidx_type_size = 0;
    int values_type_size = 0;

    MPI_Offset rowptr_nbytes = -1;
    MPI_Offset colidx_nbytes = -1;
    MPI_Offset values_nbytes = -1;

    CATCH_MPI_ERROR(MPI_Type_size(rowptr_type, &rowptr_type_size));
    CATCH_MPI_ERROR(MPI_Type_size(colidx_type, &colidx_type_size));
    CATCH_MPI_ERROR(MPI_Type_size(values_type, &values_type_size));

    // printf("grows=%ld\n", (long)crs->grows);

    {
        // Write rowptr
        CATCH_MPI_ERROR(MPI_File_open(comm, rowptr_path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));

        MPI_File_set_size(file, (crs->grows + 1) * rowptr_type_size);

        CATCH_MPI_ERROR(MPI_File_write_at_all(file,
                                              crs->rowoffset * rowptr_type_size,
                                              crs->rowptr,
                                              crs->lrows + (rank == size - 1),
                                              rowptr_type,
                                              &status));

        CATCH_MPI_ERROR(MPI_File_close(&file));
    }

    {
        // Write colidx
        CATCH_MPI_ERROR(MPI_File_open(comm, colidx_path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));

        MPI_File_set_size(file, crs->gnnz * colidx_type_size);

        CATCH_MPI_ERROR(
            MPI_File_write_at_all(file, crs->start * colidx_type_size, crs->colidx, crs->lnnz, colidx_type, &status));

        CATCH_MPI_ERROR(MPI_File_close(&file));
    }

    {
        // Write values
        CATCH_MPI_ERROR(MPI_File_open(comm, values_path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));

        MPI_File_set_size(file, crs->gnnz * values_type_size);

        CATCH_MPI_ERROR(
            MPI_File_write_at_all(file, crs->start * values_type_size, crs->values, crs->lnnz, values_type, &status));

        CATCH_MPI_ERROR(MPI_File_close(&file));
    }

    return 0;
}

int crs_write_folder(MPI_Comm comm, const char *folder, crs_t *crs) {
    assert(strlen(folder) + 11 < 1024);

    char rowptr_path[1024];
    char colidx_path[1024];
    char values_path[1024];

    sprintf(rowptr_path, "%s/rowptr.raw", folder);
    sprintf(colidx_path, "%s/colidx.raw", folder);
    sprintf(values_path, "%s/values.raw", folder);

    return crs_write(comm, rowptr_path, colidx_path, values_path, crs);
}

int crs_free(crs_t *crs) {
    free(crs->rowptr);
    free(crs->colidx);
    free(crs->values);
    crs->lrows = 0;
    crs->grows = 0;
    crs->lnnz = 0;
    crs->start = 0;
    crs->rowoffset = 0;
    return 0;
}

int crs_release(crs_t *crs) {
    crs->rowptr = 0;
    crs->colidx = 0;
    crs->values = 0;
    crs->lrows = 0;
    crs->grows = 0;
    crs->lnnz = 0;
    crs->start = 0;
    crs->rowoffset = 0;
    return 0;
}
