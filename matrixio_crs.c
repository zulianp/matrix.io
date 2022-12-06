#include "matrixio_crs.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#define CATCH_MPI_ERROR(err)                                                   \
  {                                                                            \
    if (err != MPI_SUCCESS) {                                                  \
      assert(false);                                                           \
    }                                                                          \
  }

#define MIN(a, b) ((a) < (b) ? (a) : (b))

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

int crs_read(MPI_Comm comm, const char *rowptr_path, const char *colidx_path,
             const char *values_path, MPI_Datatype rowptr_type,
             MPI_Datatype colidx_type, MPI_Datatype values_type, crs_t *crs) {
  int rank, size;

  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  MPI_File file;
  CATCH_MPI_ERROR(
      MPI_File_open(comm, rowptr_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));

  MPI_Offset rowptr_nbytes = -1;
  CATCH_MPI_ERROR(MPI_File_get_size(file, &rowptr_nbytes));

  int rowptr_type_size = 0;
  CATCH_MPI_ERROR(MPI_Type_size(rowptr_type, &rowptr_type_size));

  ptrdiff_t nrows = rowptr_nbytes / rowptr_type_size;

  if (nrows * rowptr_type_size != rowptr_nbytes) {
    if (!rank) {
      fprintf(stderr, "[Error] Wrong type specification for rowptr\n");
    }

    MPI_Abort(comm, 1);
  }

  // Remove last rowptr
  nrows -= 1;
  crs->grows = nrows;

  ptrdiff_t uniform_split = nrows / size;
  ptrdiff_t nlocal = uniform_split;
  ptrdiff_t remainder = nrows - nlocal * size;

  if (remainder > rank) {
    nlocal += 1;
  }

  ///////////////////////////////////////////////////////

  ptrdiff_t offset = rank * uniform_split;
  offset += MIN(rank, remainder);

  crs->rowptr = (void *)malloc((nlocal + 1) * rowptr_type_size);
  crs->lrows = nlocal;

  MPI_Status status;

  ///////////////////////////////////////////////////////
  // Read rowptr
  ///////////////////////////////////////////////////////

  CATCH_MPI_ERROR(MPI_File_read_at_all(file, offset * rowptr_type_size,
                                       crs->rowptr, nlocal + 1, rowptr_type,
                                       &status));
  CATCH_MPI_ERROR(MPI_File_close(&file));

  ///////////////////////////////////////////////////////
  // Read colidx
  ///////////////////////////////////////////////////////

  const char *rowptr_bytes = (const char *)crs->rowptr;

  ptrdiff_t start = to_ptrdiff_t(rowptr_type, &rowptr_bytes[0]);
  ptrdiff_t end =
      to_ptrdiff_t(rowptr_type, &(rowptr_bytes[nlocal * rowptr_type_size]));

  crs->start = start;

  ptrdiff_t nnz = end - start;
  int colidx_type_size = 0;

  crs->nnz = nnz;
  CATCH_MPI_ERROR(MPI_Type_size(colidx_type, &colidx_type_size));
  crs->colidx = (void *)malloc(nnz * colidx_type_size);

  CATCH_MPI_ERROR(
      MPI_File_open(comm, colidx_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));

  CATCH_MPI_ERROR(MPI_File_read_at_all(file, start * colidx_type_size,
                                       crs->colidx, nnz, colidx_type, &status));

  CATCH_MPI_ERROR(MPI_File_close(&file));

  ///////////////////////////////////////////////////////
  // Read values
  ///////////////////////////////////////////////////////

  int values_type_size = 0;
  CATCH_MPI_ERROR(MPI_Type_size(values_type, &values_type_size));

  crs->values = (void *)malloc(nnz * values_type_size);
  const char *values_bytes = (const char *)crs->values;

  CATCH_MPI_ERROR(
      MPI_File_open(comm, values_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));

  CATCH_MPI_ERROR(MPI_File_read_at_all(file, start * values_type_size,
                                       crs->values, nnz, values_type, &status));

  CATCH_MPI_ERROR(MPI_File_close(&file));

  ///////////////////////////////////////////////////////

  MPI_Barrier(comm);

  // for (int r = 0; r < size; r++) {
  //   if (r == rank) {
  //     printf("%d) %ld nlocal=%ld offset=%ld uniform_split=%ld remainder=%ld "
  //            "nnz=%ld, rowptr_nbytes=%d colidx_type_size=%d\n",
  //            rank, (long)nrows, nlocal, offset, uniform_split, remainder,
  //            nnz, (int)rowptr_type_size, (int)colidx_type_size);
  //     printf("start=%ld end=%ld\n", (long)start, (long)end);
  //   }

  //   MPI_Barrier(comm);
  // }

  return 0;
}

int crs_free(crs_t *crs) {
  free(crs->rowptr);
  free(crs->colidx);
  free(crs->values);
  crs->lrows = 0;
  crs->grows = 0;
  crs->nnz = 0;
  crs->start = 0;
  return 0;
}

int crs_release(crs_t *crs) {
  crs->rowptr = 0;
  crs->colidx = 0;
  crs->values = 0;
  crs->lrows = 0;
  crs->grows = 0;
  crs->nnz = 0;
  crs->start = 0;
  return 0;
}
