#include "matrixio_array.h"

#include "utils.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

int array_create_from_file(MPI_Comm comm,
                           const char *path,
                           MPI_Datatype type,
                           void **data,
                           ptrdiff_t *out_nlocal,
                           ptrdiff_t *out_nglobal) {
    // int rank, size;

    // MPI_Comm_rank(comm, &rank);
    // MPI_Comm_size(comm, &size);

    // MPI_Status status;
    // MPI_Offset nbytes;
    // MPI_File file;
    // int type_size;

    // CATCH_MPI_ERROR(MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    // CATCH_MPI_ERROR(MPI_File_get_size(file, &nbytes));
    // CATCH_MPI_ERROR(MPI_Type_size(type, &type_size));

    // ptrdiff_t n = nbytes / type_size;
    // if (n * type_size != nbytes) {
    //     assert(0);
    //     fprintf(stderr, "array_create_from_file: Wrong datatype - data pair\n");
    //     return 1;
    // }

    // ptrdiff_t uniform_split = n / size;
    // ptrdiff_t nlocal = uniform_split;
    // ptrdiff_t remainder = n - nlocal * size;

    // if (remainder > rank) {
    //     nlocal += 1;
    // }

    // *data = malloc(n * type_size);

    // ptrdiff_t offset = rank * uniform_split;
    // offset += MIN(rank, remainder);

    // CATCH_MPI_ERROR(MPI_File_read_at_all(file, offset * type_size, *data, nlocal, type, &status));
    // CATCH_MPI_ERROR(MPI_File_close(&file));

    // *out_nglobal = n;
    // *out_nlocal = nlocal;
    // return 0;

    return array_create_from_file_segmented(comm, path, type, data, INT_MAX, out_nlocal, out_nglobal);
}

int array_create_from_file_segmented(MPI_Comm comm,
                           const char *path,
                           MPI_Datatype type,
                           void **data,
                           const int segment_size,
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

    int nrounds = nlocal / segment_size;
    nrounds += nrounds * ((ptrdiff_t)segment_size) < nlocal;
    CATCH_MPI_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));

    for(int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        ptrdiff_t data_offset = i * ((ptrdiff_t)segment_size) * type_size;
        MPI_Offset byte_offset_i = (offset + i * segment_size) * type_size;

        CATCH_MPI_ERROR(MPI_File_read_at_all(file, byte_offset_i, &((char *)*data)[data_offset], segment_size_i, type, &status));
    }

    CATCH_MPI_ERROR(MPI_File_close(&file));

    *out_nglobal = n;
    *out_nlocal = nlocal;
    return 0;
}

int array_read(MPI_Comm comm, const char *path, MPI_Datatype type, void *data, ptrdiff_t nlocal, ptrdiff_t nglobal) {
    // int rank, size;

    // MPI_Comm_rank(comm, &rank);
    // MPI_Comm_size(comm, &size);

    // MPI_Status status;
    // MPI_Offset nbytes;
    // MPI_File file;
    // int type_size;

    // CATCH_MPI_ERROR(MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    // CATCH_MPI_ERROR(MPI_File_get_size(file, &nbytes));
    // CATCH_MPI_ERROR(MPI_Type_size(type, &type_size));

    // ptrdiff_t n = nbytes / type_size;
    // if (n * type_size != nbytes) {
    //     assert(0);
    //     fprintf(stderr, "array_create_from_file: Wrong datatype - data pair\n");
    //     return 1;
    // }

    // long offset = 0;
    // long nl = nlocal;

    // CATCH_MPI_ERROR(MPI_Exscan(&nl, &offset, 1, MPI_LONG, MPI_SUM, comm));

    // CATCH_MPI_ERROR(MPI_File_read_at_all(file, offset * type_size, data, nlocal, type, &status));
    // CATCH_MPI_ERROR(MPI_File_close(&file));
    // return 0;
    return array_read_segmented(comm, path, type, data, INT_MAX, nlocal, nglobal);
}

int array_read_segmented(MPI_Comm comm, const char *path, MPI_Datatype type, void *data, int segment_size, ptrdiff_t nlocal, ptrdiff_t nglobal) {
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

    int nrounds = nlocal / segment_size;
    nrounds += nrounds * ((ptrdiff_t)segment_size) < nlocal;

    CATCH_MPI_ERROR(MPI_Exscan(&nl, &offset, 1, MPI_LONG, MPI_SUM, comm));
    CATCH_MPI_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));

    for(int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        CATCH_MPI_ERROR(MPI_File_read_at_all(file, (offset + i * segment_size) * type_size, &((char *)data)[i * ((ptrdiff_t)segment_size) * type_size], segment_size_i, type, &status));
    }

    CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}

int array_write(MPI_Comm comm,
                const char *path,
                MPI_Datatype type,
                const void *data,
                ptrdiff_t nlocal,
                ptrdiff_t nglobal) {
    if(nglobal >= (ptrdiff_t)INT_MAX) {
        // Comunication free fallback by exploiting global information
        return array_write_segmented(comm, path, type, data, INT_MAX, nlocal, nglobal);
    }

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

int array_write_segmented(MPI_Comm comm,
                const char *path,
                MPI_Datatype type,
                const void *data,
                const int segment_size,
                ptrdiff_t nlocal,
                ptrdiff_t nglobal)
{
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

    int nrounds = nlocal / segment_size;
    nrounds += nrounds * ((ptrdiff_t)segment_size) < nlocal;

    if (size > 1) {
        CATCH_MPI_ERROR(MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm));
        CATCH_MPI_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));
    }

    for(int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        CATCH_MPI_ERROR(MPI_File_write_at_all(file, (offset + i * segment_size) * type_size, &((char *)data)[i * ((ptrdiff_t)segment_size) * type_size], segment_size_i, type, &status));
    }
        
    CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}


