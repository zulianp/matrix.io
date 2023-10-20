#include "matrixio_array.h"

#include "utils.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

int array_create_from_file(MPI_Comm comm,
    const char* path,
    MPI_Datatype type,
    void** data,
    ptrdiff_t* out_nlocal,
    ptrdiff_t* out_nglobal)
{
    return array_create_from_file_segmented(comm, path, type, data, INT_MAX, out_nlocal, out_nglobal);
}

int array_create_from_file_segmented(MPI_Comm comm,
    const char* path,
    MPI_Datatype type,
    void** data,
    const int segment_size,
    ptrdiff_t* out_nlocal,
    ptrdiff_t* out_nglobal)
{
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

    for (int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        ptrdiff_t data_offset = i * ((ptrdiff_t)segment_size) * type_size;
        MPI_Offset byte_offset_i = (offset + i * segment_size) * type_size;

        CATCH_MPI_ERROR(MPI_File_read_at_all(file, byte_offset_i, &((char*)*data)[data_offset], segment_size_i, type, &status));
    }

    CATCH_MPI_ERROR(MPI_File_close(&file));

    *out_nglobal = n;
    *out_nlocal = nlocal;
    return 0;
}

int array_read(MPI_Comm comm, const char* path, MPI_Datatype type, void* data, ptrdiff_t nlocal, ptrdiff_t nglobal)
{
    return array_read_segmented(comm, path, type, data, INT_MAX, nlocal, nglobal);
}

int array_read_segmented(MPI_Comm comm, const char* path, MPI_Datatype type, void* data, int segment_size, ptrdiff_t nlocal, ptrdiff_t nglobal)
{
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

    for (int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        CATCH_MPI_ERROR(MPI_File_read_at_all(file, (offset + i * segment_size) * type_size, &((char*)data)[i * ((ptrdiff_t)segment_size) * type_size], segment_size_i, type, &status));
    }

    CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}

int array_write(MPI_Comm comm,
    const char* path,
    MPI_Datatype type,
    const void* data,
    ptrdiff_t nlocal,
    ptrdiff_t nglobal)
{
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
    const char* path,
    MPI_Datatype type,
    const void* data,
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

    for (int i = 0; i < nrounds; i++) {
        int segment_size_i = MIN(segment_size, nlocal - i * ((ptrdiff_t)segment_size));
        CATCH_MPI_ERROR(MPI_File_write_at_all(file, (offset + i * segment_size) * type_size, &((char*)data)[i * ((ptrdiff_t)segment_size) * type_size], segment_size_i, type, &status));
    }

    CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}

int array_range_select(
    MPI_Comm comm,
    MPI_Datatype type,
    void* const in,
    void* const out,
    ptrdiff_t in_nlocal,
    ptrdiff_t range_start,
    ptrdiff_t range_end)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    long* parts = calloc(size + 1, sizeof(long));
    parts[rank + 1] = in_nlocal;
    CATCH_MPI_ERROR(MPI_Allreduce(MPI_IN_PLACE, &parts[1], size, MPI_LONG, MPI_SUM, comm));

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
    CATCH_MPI_ERROR(MPI_Alltoall(recv_count, 1, MPI_INT, send_count, 1, MPI_INT, comm));

    int* send_displs = malloc(size * sizeof(int));
    CATCH_MPI_ERROR(MPI_Alltoall(recv_displs, 1, MPI_INT, send_displs, 1, MPI_INT, comm));

    recv_displs[0] = 0;
    for (int r = 1; r < size; r++) {
        recv_displs[r] = recv_displs[r - 1] + recv_count[r - 1];
    }

    if (0) {
        for (int r = 0; r < size; r++) {
            MPI_Barrier(comm);

            if (r == rank) {
                if(rank == 0) {
                printf("\n-----------------\n");
                    for(int i = 0; i < size + 1; i++) {
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

    CATCH_MPI_ERROR(MPI_Alltoallv(
        in,
        send_count,
        send_displs,
        type,
        out,
        recv_count,
        recv_displs,
        type,
        comm));

    free(parts);
    free(send_count);
    free(send_displs);
    free(recv_count);
    free(recv_displs);
    return 0;
}
