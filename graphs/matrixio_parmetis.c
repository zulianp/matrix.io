#include "matrixio_crs.h"

#include "utils.h"

#include <mpi.h>
#include "parmetis.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define MPI_IDX_T MPI_INT




int crs_graph_decompose(MPI_Comm comm, crs_graph_t*crs, const int n_parts, int*const parts)
{
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    int MATRIXIO_DECOMPOSE_WEIGHTED = 0;
    MATRIXIO_READ_ENV(MATRIXIO_DECOMPOSE_WEIGHTED, atoi);

    float MATRIXIO_DECOMPOSE_WEIGHT = 1;
    MATRIXIO_READ_ENV(MATRIXIO_DECOMPOSE_WEIGHT, atof);

    idx_t *rowptr = (idx_t *)crs->rowptr;
    idx_t *colidx = (idx_t *)crs->colidx;
    const idx_t rstart = rowptr[0];
    const idx_t rowoffset = crs->rowoffset;

    idx_t *xadj = (idx_t *)malloc((crs->lrows + 1) * sizeof(idx_t));
    idx_t *adjncy = (idx_t *)malloc(crs->lnnz * sizeof(idx_t));

    xadj[0] = 0;
    for (ptrdiff_t r = 0; r < crs->lrows; r++) {
        const idx_t begin = rowptr[r] - rstart;
        const idx_t extent = rowptr[r + 1] - rowptr[r];
        const idx_t *cidx = &colidx[begin];

        xadj[r + 1] = xadj[r];
        for (idx_t k = 0; k < extent; k++) {
            if ((rowoffset + r) != cidx[k]) {
                adjncy[xadj[r + 1]] = cidx[k];
                xadj[r + 1]++;
            }
        }
    }

    idx_t *vtxdist = (idx_t *)malloc((size + 1) * sizeof(idx_t));
    MPI_Allgather(&rowoffset, 1, MPI_IDX_T, vtxdist, 1, MPI_IDX_T, comm);
    idx_t last = crs->lrows;
    MPI_Bcast(&last, 1, MPI_IDX_T, size - 1, comm);

    vtxdist[size] = vtxdist[size - 1] + last;

    if (!rank) {
        for (int r = 0; r <= size; r++) {
            printf("%d ", vtxdist[r]);
        }

        printf("\n");
    }

    assert(vtxdist[rank + 1] - vtxdist[rank] == crs->lrows);

    idx_t ncon = 1;
    idx_t *vwgt = 0;
    idx_t *adjwgt = 0;
    idx_t wgtflag = 0;
    idx_t numflag = 0;
    idx_t nparts = n_parts;

    real_t *tpwgts = (real_t *)malloc(nparts * ncon * sizeof(real_t));
    for (int r = 0; r < nparts; r++) {
        tpwgts[r] = 1. / nparts;
    }

    real_t ubvec[1] = {1.05};
    idx_t options[3] = {0, 1, 0};
    idx_t edgecut;

    int type_size;
    MPI_CATCH_ERROR(MPI_Type_size(crs->rowptr_type, &type_size));

    if (type_size != sizeof(idx_t)) {
        MPI_Abort(comm, -1);
    }

    MPI_CATCH_ERROR(MPI_Type_size(crs->colidx_type, &type_size));

    if (type_size != sizeof(idx_t)) {
        MPI_Abort(comm, -1);
    }

    if (MATRIXIO_DECOMPOSE_WEIGHTED) {
        wgtflag = 2;
        vwgt = (idx_t *)malloc(crs->lrows * sizeof(idx_t));

        for (ptrdiff_t r = 0; r < crs->lrows; r++) {
            vwgt[r] = xadj[r + 1] - xadj[r];
        }

        if (MATRIXIO_DECOMPOSE_WEIGHT != 1.f) {
            idx_t max_weight = 0;
            for (ptrdiff_t r = 0; r < crs->lrows; r++) {
                max_weight = MAX(max_weight, vwgt[r]);
            }

            for (ptrdiff_t r = 0; r < crs->lrows; r++) {
                vwgt[r] *= (MATRIXIO_DECOMPOSE_WEIGHT / max_weight);
                vwgt[r] = MAX(vwgt[r], 1);
            }
        }
    }

    int ret = ParMETIS_V3_PartKway(vtxdist,   // 0
                                   xadj,      // 1
                                   adjncy,    // 2
                                   vwgt,      // 3
                                   adjwgt,    // 4
                                   &wgtflag,  // 5
                                   &numflag,  // 6
                                   &ncon,     // 7
                                   &nparts,   // 8
                                   tpwgts,    // 9
                                   ubvec,     // 10
                                   options,   // 11
                                   &edgecut,  // 12
                                   parts,     // 13
                                   &comm);    // 14

    free(xadj);
    free(adjncy);

    if (MATRIXIO_DECOMPOSE_WEIGHTED) {
        free(vwgt);
    }

    return ret;
}

int crs_decompose(MPI_Comm comm, crs_t*crs, const int n_parts, int*const parts)
{
    crs_graph_t crs_graph;
    crs_graph_view_from_crs(crs, &crs_graph);
    int ret = crs_graph_decompose(comm, &crs_graph, n_parts, parts);
    crs_graph_release(&crs_graph);
    return ret;
}
