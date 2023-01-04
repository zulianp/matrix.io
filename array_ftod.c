#include "array_ftod.h"

void array_ftod(ptrdiff_t n, const float *const in, double *const out) {
    // if (((void *)in) == (void *)out) {
    // same memory address
    for (ptrdiff_t i = n - 1; i >= 0; --i) {
        out[i] = in[i];
    }
    // } else {
    // 	for(ptrdiff_t i = 0; i < n; ++i) {
    // 		out[i] = in[i];
    // 	}
    // }
}
