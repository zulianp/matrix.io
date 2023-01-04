#include "array_dtof.h"

void array_dtof(ptrdiff_t n, const double *const in, float *const out) {
    for (ptrdiff_t i = 0; i < n; ++i) {
        out[i] = in[i];
    }
}
