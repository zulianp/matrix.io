#include "utils.h"

#include <string.h>

ptrdiff_t to_ptrdiff_t(MPI_Datatype type, const char *data) {
    if (type == MPI_INT) {
        return *((int *)data);
    }

    if (type == MPI_LONG) {
        return *((long *)data);
    }

    if (type == MPI_UNSIGNED_LONG) {
        return *((unsigned long *)data);
    }

    // Add missing cases!!?
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 666;
}

double to_double(MPI_Datatype type, const char *data) {
    if (type == MPI_FLOAT) {
        return *((float *)data);
    }

    if (type == MPI_DOUBLE) {
        return *((double *)data);
    }

    if (type == MPI_INT) {
        return *((int *)data);
    }

    if (type == MPI_LONG) {
        return *((long *)data);
    }

    // Add missing cases!!?
    assert(0);
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 666;
}

MPI_Datatype string_to_mpi_datatype(const char *name) {
    if (!strcmp(name, "float")) {
        return MPI_FLOAT;
    }

    if (!strcmp(name, "double")) {
        return MPI_DOUBLE;
    }

    if (!strcmp(name, "int")) {
        return MPI_INT;
    }

    if (!strcmp(name, "long")) {
        return MPI_LONG;
    }

    MPI_Abort(MPI_COMM_WORLD, 1);
    return MPI_CHAR;
}
