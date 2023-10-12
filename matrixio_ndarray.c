#include "matrixio_ndarray.h"

#include "utils.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

int ndarray_read(MPI_Comm comm,
                 const char *path,
                 MPI_Datatype type,
                 int ndims,
                 void *data,
                 ptrdiff_t *const nlocal,
                 const ptrdiff_t *const nglobal) {
    return ndarray_read_segmented(comm, path, type, ndims, data, INT_MAX, nlocal, nglobal);
}

int ndarray_read_segmented(MPI_Comm comm,
                           const char *path,
                           MPI_Datatype type,
                           int ndims,
                           void *data,
                           int segment_size,
                           ptrdiff_t *const nlocal,
                           const ptrdiff_t *const nglobal) {
    // nlocal is ignored and overridden for now
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

    ptrdiff_t ntotal = nbytes / type_size;
    ptrdiff_t ntotal_user = 1;

    for (int d = 0; d < ndims; d++) {
        ntotal_user *= nglobal[d];
    }

    if (ntotal_user != ntotal) {
        assert(0);
        fprintf(stderr, "ndarray_create_from_file: Wrong input sizes!\n");
        return 1;
    }

    if (ntotal * type_size != nbytes) {
        assert(0);
        fprintf(stderr, "ndarray_create_from_file: Wrong datatype - data pair!\n");
        return 1;
    }

    ptrdiff_t nlast = nglobal[ndims - 1];
    ptrdiff_t nlast_uniform_split = nlast / size;
    ptrdiff_t nlast_local = nlast_uniform_split;
    ptrdiff_t nlast_remainder = ntotal - nlast_local * size;

    if (nlast_remainder > rank) {
        nlast_local += 1;
    }

    for (int d = 0; d < ndims - 1; d++) {
        nlocal[d] = nglobal[d];
    }

    nlocal[ndims - 1] = nlast_local;

    ptrdiff_t stride = ntotal / nlast;
    long offset = 0;
    long nl = nlast_local * stride;

    int nrounds = nl / segment_size;
    nrounds += nrounds * ((ptrdiff_t)segment_size) < nl;

    CATCH_MPI_ERROR(MPI_Exscan(&nl, &offset, 1, MPI_LONG, MPI_SUM, comm));
    CATCH_MPI_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));

    for (int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nl - i * ((ptrdiff_t)segment_size));
        CATCH_MPI_ERROR(MPI_File_read_at_all(file,
                                             (offset + i * segment_size) * type_size,
                                             &((char *)data)[i * ((ptrdiff_t)segment_size) * type_size],
                                             segment_size_i,
                                             type,
                                             &status));
    }

    CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}

int ndarray_write(MPI_Comm comm,
                  const char *path,
                  MPI_Datatype type,
                  int ndims,
                  const void *data,
                  const ptrdiff_t *const nlocal,
                  const ptrdiff_t *const nglobal) {
    ptrdiff_t ntotal = 1;
    ptrdiff_t ntotal_local = 1;

    for (int d = 0; d < ndims; d++) {
        ntotal_local *= nlocal[d];
        ntotal *= nglobal[d];
    }

    if (ntotal >= (ptrdiff_t)INT_MAX) {
        // Comunication free fallback by exploiting global information
        return ndarray_write_segmented(comm, path, type, ndims, data, INT_MAX, nlocal, nglobal);
    }

    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    assert(ntotal_local <= ntotal);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size;

    CATCH_MPI_ERROR(MPI_Type_size(type, &type_size));
    nbytes = ntotal * type_size;

    CATCH_MPI_ERROR(MPI_File_open(comm, path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));
    MPI_File_set_size(file, nbytes);

    long lnl = ntotal_local;
    long offset = 0;

    if (size > 1) {
        CATCH_MPI_ERROR(MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm));
    }

    CATCH_MPI_ERROR(MPI_File_write_at_all(file, offset * type_size, data, ntotal_local, type, &status));

    CATCH_MPI_ERROR(MPI_File_close(&file));

    return 0;
}

int ndarray_write_segmented(MPI_Comm comm,
                            const char *path,
                            MPI_Datatype type,
                            int ndims,
                            const void *data,
                            const int segment_size,
                            const ptrdiff_t *const nlocal,
                            const ptrdiff_t *const nglobal) {
    // int rank, size;

    // MPI_Comm_rank(comm, &rank);
    // MPI_Comm_size(comm, &size);

    // assert(nlocal <= nglobal);

    // MPI_Status status;
    // MPI_Offset nbytes;
    // MPI_File file;
    // int type_size;

    // CATCH_MPI_ERROR(MPI_Type_size(type, &type_size));
    // nbytes = nglobal * type_size;

    // CATCH_MPI_ERROR(MPI_File_open(comm, path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));
    // MPI_File_set_size(file, nbytes);

    // long lnl = nlocal;
    // long offset = 0;

    // int nrounds = nlocal / segment_size;
    // nrounds += nrounds * ((ptrdiff_t)segment_size) < nlocal;

    // if (size > 1) {
    //     CATCH_MPI_ERROR(MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm));
    //     CATCH_MPI_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));
    // }

    // for (int i = 0; i < nrounds; i++) {
    //     int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
    //     CATCH_MPI_ERROR(MPI_File_write_at_all(file,
    //                                           (offset + i * segment_size) * type_size,
    //                                           &((char *)data)[i * ((ptrdiff_t)segment_size) * type_size],
    //                                           segment_size_i,
    //                                           type,
    //                                           &status));
    // }

    // CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}
