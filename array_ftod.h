#ifndef ARRAY_FTOD_H
#define ARRAY_FTOD_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void array_ftod(ptrdiff_t n, const float* in, double* const out);

#ifdef __cplusplus
}
#endif

#endif
