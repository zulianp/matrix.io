#include "matrixio_array.h"

#include "utils.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int MATRIXIO_SEGMENT_CONVERT_MAX_DEFAULT = 1000000;

int array_create_from_file(MPI_Comm comm,
                           const char* path,
                           MPI_Datatype type,
                           void** data,
                           ptrdiff_t* out_nlocal,
                           ptrdiff_t* out_nglobal) {
    return array_create_from_file_segmented(comm, path, type, data, INT_MAX, out_nlocal, out_nglobal);
}

int array_create_from_file_segmented(MPI_Comm comm,
                                     const char* path,
                                     MPI_Datatype type,
                                     void** data,
                                     const int segment_size,
                                     ptrdiff_t* out_nlocal,
                                     ptrdiff_t* out_nglobal) {
    assert(!mpi_type_file_compatible(type, path));

    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size;

    MPI_CATCH_ERROR(MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    MPI_CATCH_ERROR(MPI_File_get_size(file, &nbytes));
    MPI_CATCH_ERROR(MPI_Type_size(type, &type_size));

    ptrdiff_t n = nbytes / type_size;
    if (n * type_size != nbytes) {
        assert(0);
        fprintf(stderr, "array_create_from_file: Wrong datatype - data pair\n");
        fflush(stderr);
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
    MPI_CATCH_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));

    for (int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        ptrdiff_t data_offset = i * ((ptrdiff_t)segment_size) * type_size;
        MPI_Offset byte_offset_i = (offset + i * segment_size) * type_size;

        MPI_CATCH_ERROR(
            MPI_File_read_at_all(file, byte_offset_i, &((char*)*data)[data_offset], segment_size_i, type, &status));
    }

    MPI_CATCH_ERROR(MPI_File_close(&file));

    *out_nglobal = n;
    *out_nlocal = nlocal;
    return 0;
}

int array_read(MPI_Comm comm, const char* path, MPI_Datatype type, void* data, ptrdiff_t nlocal, ptrdiff_t nglobal) {
    return array_read_segmented(comm, path, type, data, INT_MAX, nlocal, nglobal);
}

int array_read_segmented(MPI_Comm comm,
                         const char* path,
                         MPI_Datatype type,
                         void* data,
                         int segment_size,
                         ptrdiff_t nlocal,
                         ptrdiff_t nglobal) {
    assert(!mpi_type_file_compatible(type, path));

    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size;

    MPI_CATCH_ERROR(MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    MPI_CATCH_ERROR(MPI_File_get_size(file, &nbytes));
    MPI_CATCH_ERROR(MPI_Type_size(type, &type_size));

    ptrdiff_t n = nbytes / type_size;
    if (n * type_size != nbytes) {
        assert(0);
        fprintf(stderr, "array_read_segmented: Wrong datatype - data pair\n");
        return 1;
    }

    long offset = 0;
    long nl = nlocal;

    int nrounds = nlocal / segment_size;
    nrounds += nrounds * ((ptrdiff_t)segment_size) < nlocal;

    MPI_CATCH_ERROR(MPI_Exscan(&nl, &offset, 1, MPI_LONG, MPI_SUM, comm));
    MPI_CATCH_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));

    for (int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        MPI_CATCH_ERROR(MPI_File_read_at_all(file,
                                             (offset + i * segment_size) * type_size,
                                             &((char*)data)[i * ((ptrdiff_t)segment_size) * type_size],
                                             segment_size_i,
                                             type,
                                             &status));
    }

    MPI_CATCH_ERROR(MPI_File_close(&file));
    return 0;
}

int array_read_convert(MPI_Comm comm,
                       const char* path,
                       MPI_Datatype type,
                       void* data,
                       ptrdiff_t nlocal,
                       ptrdiff_t nglobal) {
    int type_size;
    MPI_CATCH_ERROR(MPI_Type_size(type, &type_size));

    MPI_Datatype file_type = mpi_type_from_file_extension(path);
    if (file_type == type || file_type == MPI_DATATYPE_NULL) {
        return array_read_segmented(comm, path, type, data, INT_MAX, nlocal, nglobal);
    }

    int file_type_size;
    MPI_CATCH_ERROR(MPI_Type_size(file_type, &file_type_size));

    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int MATRIXIO_SEGMENT_CONVERT_MAX = MATRIXIO_SEGMENT_CONVERT_MAX_DEFAULT;
    MATRIXIO_READ_ENV(MATRIXIO_SEGMENT_CONVERT_MAX, atoi);
    const int segment_size = MIN(MATRIXIO_SEGMENT_CONVERT_MAX, nlocal);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;

    MPI_CATCH_ERROR(MPI_File_open(comm, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    MPI_CATCH_ERROR(MPI_File_get_size(file, &nbytes));

    ptrdiff_t n = nbytes / file_type_size;
    if (n * file_type_size != nbytes) {
        assert(0);
        fprintf(stderr, "array_create_from_file: Wrong datatype - data pair\n");
        fflush(stderr);
        return 1;
    }

    void* buffer = 0;

    if (file_type_size > type_size) {
        printf("ALLOCATING BUFFER!\n");
        buffer = malloc(segment_size * file_type_size);
    } else {
        buffer = data;
    }

    long lnl = nlocal;
    long offset = 0;
    MPI_CATCH_ERROR(MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm));

    int nrounds = nlocal / segment_size;
    nrounds += nrounds * ((ptrdiff_t)segment_size) < nlocal;
    MPI_CATCH_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));

    for (int i = 0; i < nrounds; i++) {
        const ptrdiff_t data_offset = i * ((ptrdiff_t)segment_size) * type_size;
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));

        if (file_type_size <= type_size) {
            buffer = &((char*)data)[data_offset];
        }

        MPI_Offset byte_offset_i = (offset + i * segment_size) * file_type_size;
        MPI_CATCH_ERROR(MPI_File_read_at_all(file, byte_offset_i, buffer, segment_size_i, file_type, &status));
        array_convert(segment_size_i, file_type, buffer, type, &((char*)data)[data_offset]);
    }

    if (file_type_size > type_size) {
        free(buffer);
    }

    MPI_CATCH_ERROR(MPI_File_close(&file));
    return 0;
}

int array_write(MPI_Comm comm,
                const char* path,
                MPI_Datatype type,
                const void* data,
                ptrdiff_t nlocal,
                ptrdiff_t nglobal) {
    if (nglobal >= (ptrdiff_t)INT_MAX) {
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

    MPI_CATCH_ERROR(MPI_Type_size(type, &type_size));
    nbytes = nglobal * type_size;

    MPI_CATCH_ERROR(MPI_File_open(comm, path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));
    MPI_File_set_size(file, nbytes);

    long lnl = nlocal;
    long offset = 0;

    if (size > 1) {
        MPI_CATCH_ERROR(MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm));
    }

    MPI_CATCH_ERROR(MPI_File_write_at_all(file, offset * type_size, data, nlocal, type, &status));

    MPI_CATCH_ERROR(MPI_File_close(&file));

    return 0;
}

int array_write_segmented(MPI_Comm comm,
                          const char* path,
                          MPI_Datatype type,
                          const void* data,
                          const int segment_size,
                          ptrdiff_t nlocal,
                          ptrdiff_t nglobal) {
    assert(!mpi_type_file_compatible(type, path));
    
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    assert(nlocal <= nglobal);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size;

    MPI_CATCH_ERROR(MPI_Type_size(type, &type_size));
    nbytes = nglobal * type_size;

    MPI_CATCH_ERROR(MPI_File_open(comm, path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));
    MPI_File_set_size(file, nbytes);

    long lnl = nlocal;
    long offset = 0;

    int nrounds = nlocal / segment_size;
    nrounds += nrounds * ((ptrdiff_t)segment_size) < nlocal;

    if (size > 1) {
        MPI_CATCH_ERROR(MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm));
        MPI_CATCH_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));
    }

    for (int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        MPI_CATCH_ERROR(MPI_File_write_at_all(file,
                                              (offset + i * segment_size) * type_size,
                                              &((char*)data)[i * ((ptrdiff_t)segment_size) * type_size],
                                              segment_size_i,
                                              type,
                                              &status));
    }

    MPI_CATCH_ERROR(MPI_File_close(&file));
    return 0;
}

int array_range_select(MPI_Comm comm,
                       MPI_Datatype type,
                       void* const in,
                       void* const out,
                       ptrdiff_t in_nlocal,
                       ptrdiff_t range_start,
                       ptrdiff_t range_end) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    long* parts = calloc(size + 1, sizeof(long));
    parts[rank + 1] = in_nlocal;
    MPI_CATCH_ERROR(MPI_Allreduce(MPI_IN_PLACE, &parts[1], size, MPI_LONG, MPI_SUM, comm));

    for (int r = 1; r < size; r++) {
        parts[r + 1] += parts[r];
    }

    int* recv_count = calloc(size, sizeof(int));
    int* recv_displs = malloc(size * sizeof(int));

    for (int r = 0; r < size; r++) {
        ptrdiff_t bmin = MAX(parts[r], range_start);
        ptrdiff_t bmax = MIN(parts[r + 1], range_end);
        recv_count[r] = MAX(0, bmax - bmin);

        // From which offset process r needs to send data to this rank
        recv_displs[r] = bmin - parts[r];
    }

    int* send_count = malloc(size * sizeof(int));
    MPI_CATCH_ERROR(MPI_Alltoall(recv_count, 1, MPI_INT, send_count, 1, MPI_INT, comm));

    int* send_displs = malloc(size * sizeof(int));
    MPI_CATCH_ERROR(MPI_Alltoall(recv_displs, 1, MPI_INT, send_displs, 1, MPI_INT, comm));

    recv_displs[0] = 0;
    for (int r = 1; r < size; r++) {
        recv_displs[r] = recv_displs[r - 1] + recv_count[r - 1];
    }

    if (0) {
        for (int r = 0; r < size; r++) {
            MPI_Barrier(comm);

            if (r == rank) {
                if (rank == 0) {
                    printf("\n-----------------\n");
                    for (int i = 0; i < size + 1; i++) {
                        printf("%ld ", parts[i]);
                    }
                }

                printf("\n-----------------\n");
                printf("[%d]\n", rank);
                printf("[%ld, %ld) size = %ld\n", range_start, range_end, range_end - range_start);
                printf("\nsend_counts:\n");
                for (int i = 0; i < size; i++) {
                    printf("%d ", send_count[i]);
                }

                printf("\nsend_displs:\n");
                for (int i = 0; i < size; i++) {
                    printf("%d ", send_displs[i]);
                }

                printf("\nrecv_counts:\n");
                for (int i = 0; i < size; i++) {
                    printf("%d ", recv_count[i]);
                }

                printf("\nrecv_displs:\n");
                for (int i = 0; i < size; i++) {
                    printf("%d ", recv_displs[i]);
                }

                printf("\n-----------------\n");
            }

            fflush(stdout);
            MPI_Barrier(comm);
        }
    }

    MPI_CATCH_ERROR(MPI_Alltoallv(in, send_count, send_displs, type, out, recv_count, recv_displs, type, comm));

    free(parts);
    free(send_count);
    free(send_displs);
    free(recv_count);
    free(recv_displs);
    return 0;
}

#define ARRAY_CONVERT_(from_mpi_type_, from_type_, to_mpi_type_, to_type_) \
    do {                                                                   \
        if (from_mpi_type_ == from_type) {                                 \
            from_type_* d_from = (from_type_*)from;                        \
            if (to_mpi_type_ == to_type) {                                 \
                to_type_* d_to = (to_type_*)to;                            \
                for (ptrdiff_t i = 0; i < n; i++) {                        \
                    d_to[i] = d_from[i];                                   \
                }                                                          \
                return 0;                                                  \
            }                                                              \
        }                                                                  \
    } while (0)

// For supporting in-place conversion
#define ARRAY_CONVERT_REVERSE_(from_mpi_type_, from_type_, to_mpi_type_, to_type_) \
    do {                                                                           \
        if (from_mpi_type_ == from_type) {                                         \
            from_type_* d_from = (from_type_*)from;                                \
            if (to_mpi_type_ == to_type) {                                         \
                to_type_* d_to = (to_type_*)to;                                    \
                for (ptrdiff_t i = n - 1; i >= 0; --i) {                           \
                    d_to[i] = d_from[i];                                           \
                }                                                                  \
                return 0;                                                          \
            }                                                                      \
        }                                                                          \
    } while (0)

int array_convert(const ptrdiff_t n, MPI_Datatype from_type, const void* from, MPI_Datatype to_type, void* to) {
    int from_size;
    MPI_CATCH_ERROR(MPI_Type_size(from_type, &from_size));

    if (from_type == to_type) {
        memcpy(to, from, n * from_size);
        return 0;
    }

    // Floating points
    ARRAY_CONVERT_REVERSE_(MPI_FLOAT, float, MPI_DOUBLE, double);
    ARRAY_CONVERT_(MPI_DOUBLE, double, MPI_FLOAT, float);

    // Indices
    ARRAY_CONVERT_(MPI_LONG, long, MPI_INT, int);
    ARRAY_CONVERT_(MPI_LONG, long, MPI_INT32_T, int32_t);
    ARRAY_CONVERT_(MPI_LONG, long, MPI_INT64_T, int64_t);

    ARRAY_CONVERT_(MPI_INT, int, MPI_INT32_T, int32_t);
    ARRAY_CONVERT_REVERSE_(MPI_INT, int, MPI_LONG, long);
    ARRAY_CONVERT_REVERSE_(MPI_INT, int, MPI_INT64_T, int64_t);

    ARRAY_CONVERT_(MPI_INT32_T, int32_t, MPI_INT, int);
    ARRAY_CONVERT_REVERSE_(MPI_INT32_T, int32_t, MPI_LONG, long);
    ARRAY_CONVERT_REVERSE_(MPI_INT32_T, int32_t, MPI_INT64_T, int64_t);

    ARRAY_CONVERT_(MPI_INT64_T, int64_t, MPI_INT, int);
    ARRAY_CONVERT_(MPI_INT64_T, int64_t, MPI_LONG, long);
    ARRAY_CONVERT_(MPI_INT64_T, int64_t, MPI_INT32_T, int32_t);

    fprintf(stderr, "array_convert: Invalid convertion!\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
}

int array_write_convert(MPI_Comm comm,
                        const char* path,
                        MPI_Datatype type,
                        const void* data,
                        ptrdiff_t nlocal,
                        ptrdiff_t nglobal) {
    MPI_Datatype file_type = mpi_type_from_file_extension(path);
    if (file_type == type || file_type == MPI_DATATYPE_NULL) {
        return array_write_segmented(comm, path, type, data, INT_MAX, nlocal, nglobal);
    }

    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    assert(nlocal <= nglobal);

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size, file_type_size;

    MPI_CATCH_ERROR(MPI_Type_size(type, &type_size));
    MPI_CATCH_ERROR(MPI_Type_size(file_type, &file_type_size));

    nbytes = nglobal * file_type_size;
    MPI_CATCH_ERROR(MPI_File_open(comm, path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));
    MPI_File_set_size(file, nbytes);

    long lnl = nlocal;
    long offset = 0;

    int MATRIXIO_SEGMENT_CONVERT_MAX = MATRIXIO_SEGMENT_CONVERT_MAX_DEFAULT;
    MATRIXIO_READ_ENV(MATRIXIO_SEGMENT_CONVERT_MAX, atoi);
    const ptrdiff_t segment_size = MIN(MATRIXIO_SEGMENT_CONVERT_MAX, nlocal);
    int nrounds = nlocal / segment_size;
    nrounds += nrounds * segment_size < nlocal;

    if (size > 1) {
        MPI_CATCH_ERROR(MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm));
        MPI_CATCH_ERROR(MPI_Exscan(MPI_IN_PLACE, &nrounds, 1, MPI_INT, MPI_MAX, comm));
    }

    void* buffer = malloc(segment_size * file_type_size);
    for (int i = 0; i < nrounds; i++) {
        const ptrdiff_t segment_size_i = MIN(segment_size, nlocal - i * segment_size);

        array_convert(
            segment_size_i, type, &((char*)data)[i * segment_size * type_size], file_type, buffer);

        MPI_CATCH_ERROR(MPI_File_write_at_all(
            file, (offset + i * segment_size) * file_type_size, buffer, segment_size_i, file_type, &status));
    }

    free(buffer);
    MPI_CATCH_ERROR(MPI_File_close(&file));

    return 0;
}

#undef ARRAY_CONVERT_
