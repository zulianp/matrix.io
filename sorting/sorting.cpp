#include "sorting.h"

#include <mpi.h>
#include <stdio.h>
#include <algorithm>
#include <cstring>

template <typename KeyType>
inline static void local_argsort_impl(ptrdiff_t size, const KeyType *keys, ptrdiff_t *idx) {
    for (ptrdiff_t i = 0; i < size; i++) {
        idx[i] = i;
    }

    auto comp = [&keys](const ptrdiff_t l, const ptrdiff_t r) -> bool { return keys[l] < keys[r]; };
    std::sort(idx, idx + size, comp);
}

extern "C" void local_argsort(ptrdiff_t size, MPI_Datatype key_type, void *const keys, ptrdiff_t *idx) {
    if (key_type == MPI_INT) {
        local_argsort_impl<int>(size, (int *)keys, idx);
    } else if (key_type == MPI_LONG) {
        local_argsort_impl<long>(size, (long *)keys, idx);
    } else {
        fprintf(stderr, "Unsupported key type!\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

extern "C" void gather_remap(ptrdiff_t size, int type_size, const ptrdiff_t *const gather_idx, void *const input, void *const output)
{
    for(ptrdiff_t i =0; i < size; i++) {
        memcpy(&((char *)output)[i*type_size], &((char *)input)[gather_idx[i]*type_size], type_size);
    }
}