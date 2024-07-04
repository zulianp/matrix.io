#include "matrixio_crs.h"

#include "utils.h"

#include <assert.h>
#include <glob.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PATH_LENGTH 4096

int crs_alloc_same(const crs_t *const tpl, crs_t *result) {
    result->grows = tpl->grows;
    result->lrows = tpl->lrows;
    result->lnnz = tpl->lnnz;
    result->gnnz = tpl->gnnz;
    result->start = tpl->start;
    result->rowoffset = tpl->rowoffset;

    result->rowptr_type = tpl->rowptr_type;
    result->colidx_type = tpl->colidx_type;
    result->values_type = tpl->values_type;

    int rowptr_type_size = 0;
    CATCH_MPI_ERROR(MPI_Type_size(tpl->rowptr_type, &rowptr_type_size));

    int colidx_type_size = 0;
    CATCH_MPI_ERROR(MPI_Type_size(tpl->colidx_type, &colidx_type_size));

    int values_type_size = 0;
    CATCH_MPI_ERROR(MPI_Type_size(tpl->values_type, &values_type_size));

    result->rowptr = malloc(tpl->lrows * rowptr_type_size);
    result->colidx = malloc(tpl->lnnz * colidx_type_size);
    result->values = malloc(tpl->lnnz * values_type_size);
    return 0;
}

int crs_read_str(MPI_Comm comm,
                 const char *rowptr_path,
                 const char *colidx_path,
                 const char *values_path,
                 const char *rowptr_type_str,
                 const char *colidx_type_str,
                 const char *values_type_str,
                 crs_t *crs) {
    MPI_Datatype rowptr_type = string_to_mpi_datatype(rowptr_type_str);
    MPI_Datatype colidx_type = string_to_mpi_datatype(colidx_type_str);
    MPI_Datatype values_type = string_to_mpi_datatype(values_type_str);
    return crs_read(comm, rowptr_path, colidx_path, values_path, rowptr_type, colidx_type, values_type, crs);
}

int crs_graph_read_AoS_block(MPI_Comm comm,
                             const char *rowptr_path,
                             const char *colidx_path,
                             MPI_Datatype rowptr_type,
                             MPI_Datatype colidx_type,
                             int block_size,
                             crs_graph_t *crs) {
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_File file;
    MPI_Offset rowptr_nbytes = -1;
    int rowptr_type_size = 0;

    CATCH_MPI_ERROR(MPI_File_open(comm, rowptr_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));
    CATCH_MPI_ERROR(MPI_File_get_size(file, &rowptr_nbytes));
    CATCH_MPI_ERROR(MPI_Type_size(rowptr_type, &rowptr_type_size));

    ptrdiff_t nrows = rowptr_nbytes / rowptr_type_size;

    if (nrows * rowptr_type_size != rowptr_nbytes) {
        if (!rank) {
            fprintf(stderr,
                    "[Error] Wrong type specification for rowptr (%ld = %ld / %d) in file %s\n",
                    (long)nrows,
                    (long)rowptr_nbytes,
                    (int)rowptr_type_size,
                    rowptr_path);

            fflush(stderr);
        }

        MPI_Abort(comm, 1);
    }

    // Remove last rowptr
    nrows -= 1;
    crs->grows = nrows;

    ptrdiff_t uniform_split = block_size * (nrows / (size * block_size));
    ptrdiff_t nlocal = uniform_split;
    ptrdiff_t remainder = nrows - nlocal * size;

    if (remainder > rank * block_size) {
        nlocal += block_size;
    }

#ifndef NDEBUG
    long ntotal = nlocal;
    MPI_Allreduce(MPI_IN_PLACE, &ntotal, 1, MPI_LONG, MPI_SUM, comm);
    if(ntotal != nrows) {
        printf("ntotal != nrows, %ld != %ld!\n", ntotal, nrows);    
        fflush(stdout);
    }
    assert(ntotal == nrows);
    
#endif

    ///////////////////////////////////////////////////////

    ptrdiff_t offset = rank * uniform_split;
    offset += MIN(rank * block_size, remainder);

    char *rowptr = (char *)malloc((nlocal + 1) * rowptr_type_size);

    MPI_Status status;

    ///////////////////////////////////////////////////////
    // Read rowptr
    ///////////////////////////////////////////////////////

    CATCH_MPI_ERROR(MPI_File_read_at_all(file, offset * rowptr_type_size, rowptr, nlocal + 1, rowptr_type, &status));
    CATCH_MPI_ERROR(MPI_File_close(&file));

    ///////////////////////////////////////////////////////
    // Read colidx
    ///////////////////////////////////////////////////////

    MPI_Offset gnnz_bytes = -1;
    ptrdiff_t start = to_ptrdiff_t(rowptr_type, &rowptr[0]);
    ptrdiff_t end = to_ptrdiff_t(rowptr_type, &(rowptr[nlocal * rowptr_type_size]));

    crs->start = start;

    ptrdiff_t nnz = end - start;
    int colidx_type_size = 0;

    crs->lnnz = nnz;
    CATCH_MPI_ERROR(MPI_Type_size(colidx_type, &colidx_type_size));

    char *colidx = (char *)malloc(nnz * colidx_type_size);

    CATCH_MPI_ERROR(MPI_File_open(comm, colidx_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));

    CATCH_MPI_ERROR(MPI_File_get_size(file, &gnnz_bytes));

    CATCH_MPI_ERROR(MPI_File_read_at_all(file, start * colidx_type_size, colidx, nnz, colidx_type, &status));

    CATCH_MPI_ERROR(MPI_File_close(&file));

    ///////////////////////////////////////////////////////

    crs->rowptr = rowptr;
    crs->colidx = colidx;

    crs->grows = nrows;
    crs->lrows = nlocal;
    crs->gnnz = gnnz_bytes / colidx_type_size;
    crs->lnnz = nnz;
    crs->start = start;

    crs->rowptr_type = rowptr_type;
    crs->colidx_type = colidx_type;
    crs->rowoffset = offset;

    return 0;
}

int crs_graph_read(MPI_Comm comm,
                   const char *rowptr_path,
                   const char *colidx_path,
                   MPI_Datatype rowptr_type,
                   MPI_Datatype colidx_type,
                   crs_graph_t *crs) {
    return crs_graph_read_AoS_block(comm, rowptr_path, colidx_path, rowptr_type, colidx_type, 1, crs);
}

int crs_graph_read_values(MPI_Comm comm,
                          const crs_graph_t *const crs,
                          const char *values_path,
                          MPI_Datatype values_type,
                          matrixio_byte_t *const values) {
    ///////////////////////////////////////////////////////
    // Read values
    ///////////////////////////////////////////////////////

    MPI_File file;
    MPI_Status status;

    int values_type_size = 0;
    CATCH_MPI_ERROR(MPI_Type_size(values_type, &values_type_size));

    CATCH_MPI_ERROR(MPI_File_open(comm, values_path, MPI_MODE_RDONLY, MPI_INFO_NULL, &file));

    CATCH_MPI_ERROR(MPI_File_read_at_all(file, crs->start * values_type_size, values, crs->lnnz, values_type, &status));

    CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}

int crs_graph_to_crs(crs_graph_t *const graph, crs_t *const crs) {
    crs->rowptr = graph->rowptr;
    crs->colidx = graph->colidx;

    crs->grows = graph->grows;
    crs->lrows = graph->lrows;
    crs->lnnz = graph->lnnz;
    crs->gnnz = graph->gnnz;
    crs->start = graph->start;
    crs->rowoffset = graph->rowoffset;

    crs->rowptr_type = graph->rowptr_type;
    crs->colidx_type = graph->colidx_type;

    // Make graph invalid
    crs_graph_release(graph);
    return 0;
}

int crs_graph_view_from_crs(crs_t *const crs, crs_graph_t *const graph) {
    graph->rowptr = crs->rowptr;
    graph->colidx = crs->colidx;

    graph->grows = crs->grows;
    graph->lrows = crs->lrows;
    graph->lnnz = crs->lnnz;
    graph->gnnz = crs->gnnz;
    graph->start = crs->start;
    graph->rowoffset = crs->rowoffset;

    graph->rowptr_type = crs->rowptr_type;
    graph->colidx_type = crs->colidx_type;
    return 0;
}

int crs_graph_view_from_block_crs(block_crs_t *const crs, crs_graph_t *const graph) {
    graph->rowptr = crs->rowptr;
    graph->colidx = crs->colidx;

    graph->grows = crs->grows;
    graph->lrows = crs->lrows;
    graph->lnnz = crs->lnnz;
    graph->gnnz = crs->gnnz;
    graph->start = crs->start;
    graph->rowoffset = crs->rowoffset;

    graph->rowptr_type = crs->rowptr_type;
    graph->colidx_type = crs->colidx_type;
    return 0;
}

int crs_graph_to_block_crs(crs_graph_t *const graph, block_crs_t *const crs) {
    crs->rowptr = graph->rowptr;
    crs->colidx = graph->colidx;

    crs->grows = graph->grows;
    crs->lrows = graph->lrows;
    crs->lnnz = graph->lnnz;
    crs->gnnz = graph->gnnz;
    crs->start = graph->start;
    crs->rowoffset = graph->rowoffset;

    crs->rowptr_type = graph->rowptr_type;
    crs->colidx_type = graph->colidx_type;

    // Make graph invalid
    crs_graph_release(graph);
    return 0;
}

int crs_read(MPI_Comm comm,
             const char *rowptr_path,
             const char *colidx_path,
             const char *values_path,
             MPI_Datatype rowptr_type,
             MPI_Datatype colidx_type,
             MPI_Datatype values_type,
             crs_t *crs) {
    // crs_graph_t graph;
    // crs_graph_read(comm, rowptr_path, colidx_path, rowptr_type, colidx_type, &graph);

    // int values_type_size = 0;
    // CATCH_MPI_ERROR(MPI_Type_size(values_type, &values_type_size));
    // char *values = (char *)malloc(graph.lnnz * values_type_size);

    // crs_graph_read_values(comm, &graph, values_path, values_type, values);
    // crs_graph_to_crs(&graph, crs);

    // crs->values = values;
    // crs->values_type = values_type;
    // return 0;

    int MATRIXIO_CRS_READ_BLOCK_SIZE = 1;
    MATRIXIO_READ_ENV(MATRIXIO_CRS_READ_BLOCK_SIZE, atoi);

    if(MATRIXIO_CRS_READ_BLOCK_SIZE != 1) {
        int rank;
        MPI_Comm_rank(comm, &rank);
        if(!rank) {
            printf("MATRIXIO_CRS_READ_BLOCK_SIZE=%d\n", MATRIXIO_CRS_READ_BLOCK_SIZE);
            fflush(stdout);
        }
    }

    return crs_read_AoS_block(comm,
                              rowptr_path,
                              colidx_path,
                              values_path,
                              rowptr_type,
                              colidx_type,
                              values_type,
                              MATRIXIO_CRS_READ_BLOCK_SIZE,
                              crs);
}

int crs_read_AoS_block(MPI_Comm comm,
                       const char *rowptr_path,
                       const char *colidx_path,
                       const char *values_path,
                       MPI_Datatype rowptr_type,
                       MPI_Datatype colidx_type,
                       MPI_Datatype values_type,
                       const int block_size,
                       crs_t *crs) {
    crs_graph_t graph;
    crs_graph_read_AoS_block(comm, rowptr_path, colidx_path, rowptr_type, colidx_type, block_size, &graph);

    int values_type_size = 0;
    CATCH_MPI_ERROR(MPI_Type_size(values_type, &values_type_size));
    char *values = (char *)malloc(graph.lnnz * values_type_size);

    crs_graph_read_values(comm, &graph, values_path, values_type, values);
    crs_graph_to_crs(&graph, crs);

    crs->values = values;
    crs->values_type = values_type;
    return 0;
}

int block_crs_read(MPI_Comm comm,
                   const char *rowptr_path,
                   const char *colidx_path,
                   const char *values_pattern,
                   MPI_Datatype rowptr_type,
                   MPI_Datatype colidx_type,
                   MPI_Datatype values_type,
                   block_crs_t *crs) {
    crs_graph_t graph;
    crs_graph_read(comm, rowptr_path, colidx_path, rowptr_type, colidx_type, &graph);
    int values_type_size = 0;
    CATCH_MPI_ERROR(MPI_Type_size(values_type, &values_type_size));

    int block_size = 0;
    matrixio_byte_t **values;
    {
        glob_t gl;
        glob(values_pattern, GLOB_MARK, NULL, &gl);
        block_size = gl.gl_pathc;

        printf("block_size (%d):\n", block_size);
        for (int k = 0; k < block_size; k++) {
            printf("%s\n", gl.gl_pathv[k]);
        }

        values = (matrixio_byte_t **)malloc(block_size * sizeof(matrixio_byte_t *));

        for (int k = 0; k < block_size; k++) {
            matrixio_byte_t *v = (matrixio_byte_t *)malloc(graph.lnnz * values_type_size);
            crs_graph_read_values(comm, &graph, gl.gl_pathv[k], values_type, v);
            values[k] = v;
        }

        globfree(&gl);
    }

    crs_graph_to_block_crs(&graph, crs);

    crs->values = values;
    crs->values_type = values_type;
    crs->block_size = block_size;
    return 0;
}

int crs_read_folder(MPI_Comm comm,
                    const char *folder,
                    MPI_Datatype rowptr_type,
                    MPI_Datatype colidx_type,
                    MPI_Datatype values_type,
                    crs_t *crs) {
    assert(strlen(folder) + 11 < MAX_PATH_LENGTH);

    char rowptr_path[MAX_PATH_LENGTH];
    char colidx_path[MAX_PATH_LENGTH];
    char values_path[MAX_PATH_LENGTH];

    sprintf(rowptr_path, "%s/rowptr.raw", folder);
    sprintf(colidx_path, "%s/colidx.raw", folder);
    sprintf(values_path, "%s/values.raw", folder);

    return crs_read(comm, rowptr_path, colidx_path, values_path, rowptr_type, colidx_type, values_type, crs);
}

int crs_graph_write(MPI_Comm comm, const char *rowptr_path, const char *colidx_path, crs_graph_t *crs) {
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_File file;
    MPI_Status status;

    MPI_Datatype rowptr_type = crs->rowptr_type;
    MPI_Datatype colidx_type = crs->colidx_type;

    int rowptr_type_size = 0;
    int colidx_type_size = 0;

    CATCH_MPI_ERROR(MPI_Type_size(rowptr_type, &rowptr_type_size));
    CATCH_MPI_ERROR(MPI_Type_size(colidx_type, &colidx_type_size));

    {
        // Write rowptr
        CATCH_MPI_ERROR(MPI_File_open(comm, rowptr_path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));

        MPI_File_set_size(file, (crs->grows + 1) * rowptr_type_size);

        CATCH_MPI_ERROR(MPI_File_write_at_all(file,
                                              crs->rowoffset * rowptr_type_size,
                                              crs->rowptr,
                                              crs->lrows + (rank == size - 1),
                                              rowptr_type,
                                              &status));

        CATCH_MPI_ERROR(MPI_File_close(&file));
    }

    {
        // Write colidx
        CATCH_MPI_ERROR(MPI_File_open(comm, colidx_path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));

        MPI_File_set_size(file, crs->gnnz * colidx_type_size);

        CATCH_MPI_ERROR(
            MPI_File_write_at_all(file, crs->start * colidx_type_size, crs->colidx, crs->lnnz, colidx_type, &status));

        CATCH_MPI_ERROR(MPI_File_close(&file));
    }

    return 0;
}

int crs_graph_write_values(MPI_Comm comm,
                           const crs_graph_t *const crs,
                           const char *values_path,
                           MPI_Datatype values_type,
                           matrixio_byte_t *const values) {
    int rank, size;

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    MPI_File file;
    MPI_Status status;

    int values_type_size = 0;
    CATCH_MPI_ERROR(MPI_Type_size(values_type, &values_type_size));

    CATCH_MPI_ERROR(MPI_File_open(comm, values_path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file));

    MPI_File_set_size(file, crs->gnnz * values_type_size);

    CATCH_MPI_ERROR(
        MPI_File_write_at_all(file, crs->start * values_type_size, values, crs->lnnz, values_type, &status));

    CATCH_MPI_ERROR(MPI_File_close(&file));
    return 0;
}

int crs_write(MPI_Comm comm, const char *rowptr_path, const char *colidx_path, const char *values_path, crs_t *crs) {
    crs_graph_t graph;
    crs_graph_view_from_crs(crs, &graph);
    crs_graph_write(comm, rowptr_path, colidx_path, &graph);
    crs_graph_write_values(comm, &graph, values_path, crs->values_type, crs->values);
    crs_graph_release(&graph);
    return 0;
}

int block_crs_write(MPI_Comm comm,
                    const char *rowptr_path,
                    const char *colidx_path,
                    const char *values_format,
                    block_crs_t *crs) {
    crs_graph_t graph;
    crs_graph_view_from_block_crs(crs, &graph);
    crs_graph_write(comm, rowptr_path, colidx_path, &graph);

    char path[MAX_PATH_LENGTH];
    for (int b = 0; b < crs->block_size; b++) {
        sprintf(path, values_format, b);
        crs_graph_write_values(comm, &graph, path, crs->values_type, crs->values[b]);
    }

    crs_graph_release(&graph);
    return 0;
}

int crs_write_folder(MPI_Comm comm, const char *folder, crs_t *crs) {
    assert(strlen(folder) + 11 < MAX_PATH_LENGTH);

    char rowptr_path[MAX_PATH_LENGTH];
    char colidx_path[MAX_PATH_LENGTH];
    char values_path[MAX_PATH_LENGTH];

    sprintf(rowptr_path, "%s/rowptr.raw", folder);
    sprintf(colidx_path, "%s/colidx.raw", folder);
    sprintf(values_path, "%s/values.raw", folder);

    return crs_write(comm, rowptr_path, colidx_path, values_path, crs);
}

int crs_free(crs_t *const crs) {
    free(crs->rowptr);
    free(crs->colidx);
    free(crs->values);
    crs->lrows = 0;
    crs->grows = 0;
    crs->lnnz = 0;
    crs->start = 0;
    crs->rowoffset = 0;
    return 0;
}

int block_crs_free(block_crs_t *const crs) {
    free(crs->rowptr);
    free(crs->colidx);

    for (int b = 0; b < crs->block_size; b++) {
        free(crs->values[b]);
    }

    free(crs->values);
    crs->lrows = 0;
    crs->grows = 0;
    crs->lnnz = 0;
    crs->start = 0;
    crs->rowoffset = 0;
    crs->block_size = 0;
    return 0;
}

int crs_release(crs_t *const crs) {
    crs->rowptr = 0;
    crs->colidx = 0;
    crs->values = 0;
    crs->lrows = 0;
    crs->grows = 0;
    crs->lnnz = 0;
    crs->start = 0;
    crs->rowoffset = 0;
    return 0;
}

int crs_graph_free(crs_graph_t *const crs) {
    free(crs->rowptr);
    free(crs->colidx);

    crs->lrows = 0;
    crs->grows = 0;
    crs->lnnz = 0;
    crs->start = 0;
    crs->rowoffset = 0;
    return 0;
}

int crs_graph_release(crs_graph_t *const crs) {
    crs->rowptr = 0;
    crs->colidx = 0;

    crs->lrows = 0;
    crs->grows = 0;
    crs->lnnz = 0;
    crs->start = 0;
    crs->rowoffset = 0;
    return 0;
}

#undef MAX_PATH_LENGTH
