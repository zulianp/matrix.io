#ifndef MATRIX_IO_UTILS_H
#define MATRIX_IO_UTILS_H

#include <assert.h>
#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CATCH_MPI_ERROR(err)                                 \
    {                                                        \
        if (err != MPI_SUCCESS) {                            \
            char string_buff[4096];                          \
            int resultlen = 4096;                            \
            MPI_Error_string(err, string_buff, &resultlen);  \
            fprintf(stderr, "MPI error: %s\n", string_buff); \
            fflush(stderr);                                   \
            assert(0);                                       \
            MPI_Abort(MPI_COMM_WORLD, -1);                   \
        }                                                    \
    }

#define MATRIXIO_READ_ENV(name, conversion) \
    do {                                    \
        char *var = getenv(#name);          \
        if (var) {                          \
            name = conversion(var);         \
        }                                   \
    } while (0)

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

ptrdiff_t to_ptrdiff_t(MPI_Datatype type, const char *data);
double to_double(MPI_Datatype type, const char *data);

MPI_Datatype string_to_mpi_datatype(const char *name);
MPI_Datatype mpi_type_from_file_extension(const char *path);

#ifdef __cplusplus
}
#endif

#endif  // MATRIX_IO_UTILS_H
