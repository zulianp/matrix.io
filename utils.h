#ifndef MATRIX_IO_UTILS_H
#define MATRIX_IO_UTILS_H

#include <assert.h>
#include <mpi.h>
#include <stddef.h>
#include <stdlib.h>

#define CATCH_MPI_ERROR(err)      \
    {                             \
        if (err != MPI_SUCCESS) { \
            assert(0);            \
        }                         \
    }

#define MATRIXIO_READ_ENV(name, conversion) \
    do {                               \
        char *var = getenv(#name);     \
        if (var) {                     \
            name = conversion(var);    \
        }                              \
    } while (0)

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

ptrdiff_t to_ptrdiff_t(MPI_Datatype type, const char *data);
double to_double(MPI_Datatype type, const char *data);

MPI_Datatype string_to_mpi_datatype(const char *name);

#endif  // MATRIX_IO_UTILS_H
