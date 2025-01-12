#ifndef MATRIXIO_BASE_H
#define MATRIXIO_BASE_H

typedef char matrixio_byte_t;
#define MATRIXIO_MPI_BYTE_T MPI_CHAR

#ifdef _WIN64
#define MATRIXIO_WINDOWS 64
#else
#ifdef _WIN32
#define MATRIXIO_WINDOWS 32
#endif
#endif

#endif //MATRIXIO_BASE_H
