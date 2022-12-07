#ifndef MATRIX_IO_UTILS_H
#define MATRIX_IO_UTILS_H

#include <mpi.h>
#include <stddef.h>
#include <assert.h>

#define CATCH_MPI_ERROR(err)      \
    {                             \
        if (err != MPI_SUCCESS) { \
            assert(0);        \
        }                         \
    }

#define MIN(a, b) ((a) < (b) ? (a) : (b))

ptrdiff_t to_ptrdiff_t(MPI_Datatype type, const char *data);
double to_double(MPI_Datatype type, const char *data);

MPI_Datatype string_to_mpi_datatype(const char *name);

#endif //MATRIX_IO_UTILS_H
