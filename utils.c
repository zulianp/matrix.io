#include "utils.h"

#include <stdio.h>
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

static inline int extcmp(const char *path, const int len, const char *ext) {
    const int ext_len = strlen(ext);
    if (len < ext_len + 1) {
        return 1;
    }

    if (strcmp(&path[len - 1 - ext_len], ".float32") == 0) {
        return 0;
    }

    return 1;
}

MPI_Datatype mpi_type_from_file_extension(const char *path) {
    //
    ptrdiff_t len = strlen(path);

    // .raw (if extension is e.g., .float32.raw)
    if (extcmp(path, len, ".raw") == 0) {
        // skip
        len -= 4;
    }

    // .bin (if extension is e.g., .float32.bin)
    if (extcmp(path, len, ".bin") == 0) {
        // skip
        len -= 4;
    }

    if (extcmp(path, len, ".float32") == 0 || extcmp(path, len, ".float") == 0) {
        return MPI_FLOAT;
    }

    if (extcmp(path, len, ".float64") == 0 || extcmp(path, len, ".double") == 0) {
        return MPI_DOUBLE;
    }

    if (extcmp(path, len, ".int32") == 0) {
        return MPI_INT32_T;
    }

    if (extcmp(path, len, ".uint32") == 0) {
        return MPI_UINT32_T;
    }

    if (extcmp(path, len, ".int") == 0) {
        return MPI_INT;
    }

    if (extcmp(path, len, ".int64") == 0) {
        return MPI_INT64_T;
    }

    if (extcmp(path, len, ".long") == 0) {
        return MPI_LONG;
    }

    if (extcmp(path, len, ".ulong") == 0) {
        return MPI_UNSIGNED_LONG;
    }

    return MPI_DATATYPE_NULL;
}

int mpi_type_file_compatible(const MPI_Datatype type, const char *path) {
    MPI_Datatype file_type = mpi_type_from_file_extension(path);

    int type_size;
    MPI_CATCH_ERROR(MPI_Type_size(type, &type_size));

    int file_type_size = 0;

    if(file_type != MPI_DATATYPE_NULL) {
        MPI_CATCH_ERROR(MPI_Type_size(file_type, &file_type_size));
    }

    if (file_type == MPI_DATATYPE_NULL || file_type_size == type_size) {
        return 0;
    } else {
        fprintf(stderr, "%s non-matching extension with MPI_Datatype (%d != %d)\n", path, file_type_size, type_size);
        assert(0);
        return 1;
    }
}
