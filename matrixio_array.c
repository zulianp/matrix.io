#include "matrixio_array.h"

#include "utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

int array_create_from_file(MPI_Comm comm,
                           const char *path,
                           MPI_Datatype type,
                           void **data,
                           ptrdiff_t *out_nlocal,
                           ptrdiff_t *out_nglobal) {
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size;

    CATCH_MPI_ERROR(MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    CATCH_MPI_ERROR(MPI_File_get_size(file, &nbytes));
    CATCH_MPI_ERROR(MPI_Type_size(type, &type_size));

    ptrdiff_t n = nbytes / type_size;
    if (n * type_size != nbytes) {
        assert(0);
        fprintf(stderr, "array_create_from_file: Wrong datatype - data pair\n");
        return 1;
    }

    ptrdiff_t uniform_split = n / size;
    ptrdiff_t nlocal = uniform_split;
    ptrdiff_t remainder = n - nlocal * size;

    if (remainder > rank) {
        nlocal += 1;
    }

    *data = malloc(n * type_size);

    ptrdiff_t offset = rank * uniform_split;
    offset += MIN(rank, remainder);

    CATCH_MPI_ERROR(MPI_File_read_at_all(file, offset * type_size, *data, nlocal, type, &status));
    CATCH_MPI_ERROR(MPI_File_close(&file));

    *out_nglobal = n;
    *out_nlocal = nlocal;
    return 0;
}

int array_read(MPI_Comm comm, const char *path, MPI_Datatype type, void *data, ptrdiff_t nlocal, ptrdiff_t nglobal) {
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size;

    CATCH_MPI_ERROR(MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    CATCH_MPI_ERROR(MPI_File_get_size(file, &nbytes));
    CATCH_MPI_ERROR(MPI_Type_size(type, &type_size));

    ptrdiff_t n = nbytes / type_size;
    if (n * type_size != nbytes) {
        assert(0);
        fprintf(stderr, "array_create_from_file: Wrong datatype - data pair\n");
        return 1;
    }

    long offset = 0;
    long nl = nlocal;

    CATCH_MPI_ERROR(MPI_Exscan(&nl, &offset, 1, MPI_LONG, MPI_SUM, comm));

    CATCH_MPI_ERROR(MPI_File_read_at_all(file, offset * type_size, data, nlocal, type, &status));
    CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}

int array_write(MPI_Comm comm,
                const char *path,
                MPI_Datatype type,
                const void *data,
                ptrdiff_t nlocal,
                ptrdiff_t nglobal) {
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    assert(nlocal <= nglobal);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size;

    CATCH_MPI_ERROR(MPI_Type_size(type, &type_size));
    nbytes = nglobal * type_size;

    CATCH_MPI_ERROR(MPI_File_open(comm, path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));
    MPI_File_set_size(file, nbytes);

    long lnl = nlocal;
    long offset = 0;

    if (size > 1) {
        CATCH_MPI_ERROR(MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm));
    }

    CATCH_MPI_ERROR(MPI_File_write_at_all(file, offset * type_size, data, nlocal, type, &status));

    CATCH_MPI_ERROR(MPI_File_close(&file));

    return 0;
}
