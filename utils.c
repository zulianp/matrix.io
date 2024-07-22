#include "utils.h"

#include <string.h>

ptrdiff_t to_ptrdiff_t(MPI_Datatype type, const char *data) {
    if (type == MPI_INT) {
        return *((int *)data);
    }

    if (type == MPI_LONG) {
        return *((long *)data);
    }

    if (type == MPI_INT32_T) {
        return *((int32_t *)data);
    }

    if (type == MPI_INT64_T) {
        return *((int64_t *)data);
    }

    if (type == MPI_UNSIGNED_LONG) {
        return *((unsigned long *)data);
    }

    if (type == MPI_INT) {
        return *((int *)data);
    }

    if (type == MPI_UINT32_T) {
        return *((int *)data);
    }

    if (type == MPI_LONG) {
        return *((long *)data);
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

    if (!strcmp(name, "float32")) {
        return MPI_FLOAT;
    }

    if (!strcmp(name, "float64")) {
        return MPI_DOUBLE;
    }

    if (!strcmp(name, "int32")) {
        return MPI_INT32_T;
    }

    if (!strcmp(name, "uint32")) {
        return MPI_UINT32_T;
    }

    if (!strcmp(name, "int64")) {
        return MPI_INT64_T;
    }

    if (!strcmp(name, "ulong")) {
        return MPI_UNSIGNED_LONG;
    }

    MPI_Abort(MPI_COMM_WORLD, 1);
    return MPI_CHAR;
}

MPI_Datatype mpi_type_from_file_extension(const char *path) {
    //
    ptrdiff_t len = strlen(path);

    // .raw (if extension is e.g., .float32.raw)
    if (strcmp(&path[len - 1 - 4], ".raw") == 0) {
        // skip
        len -= 4;
    }

    // .bin (if extension is e.g., .float32.bin)
    if (strcmp(&path[len - 1 - 4], ".bin") == 0) {
        // skip
        len -= 4;
    }

    if (strcmp(&path[len - 1 - 8], ".float32") == 0 || strcmp(&path[len - 1 - 6], ".float")) {
        return MPI_FLOAT;
    }

    if (strcmp(&path[len - 1 - 8], ".float64") == 0 || strcmp(&path[len - 1 - 7], ".double") == 0) {
        return MPI_DOUBLE;
    }

    if (strcmp(&path[len - 1 - 6], ".int32") == 0) {
        return MPI_INT32_T;
    }

    if (strcmp(&path[len - 1 - 7], ".uint32") == 0) {
        return MPI_UINT32_T;
    }

    if (strcmp(&path[len - 1 - 6], ".int") == 0) {
        return MPI_INT;
    }

    if (strcmp(&path[len - 1 - 6], ".int64") == 0) {
        return MPI_INT64_T;
    }

    if (strcmp(&path[len - 1 - 5], ".long") == 0) {
        return MPI_LONG;
    }

    if (strcmp(&path[len - 1 - 6], ".ulong") == 0) {
        return MPI_UNSIGNED_LONG;
    }

    return MPI_DATATYPE_NULL;
}
