#include "matrixio_crs.h"

#include "utils.h"

#include <mpi.h>
#include "parmetis.h"

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define MPI_IDX_T MPI_INT

int decompose(MPI_Comm comm, crs_t *crs, const int n_parts, int *const parts) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    idx_t *vtxdist = (idx_t *)malloc((size + 1) * sizeof(idx_t));
    // memset(vtxdist, 0, (size + 1) * sizeof(idx_t));

    idx_t rowoffset = crs->rowoffset;
    MPI_Allgather(&rowoffset, 1, MPI_IDX_T, vtxdist, 1, MPI_IDX_T, comm);

    idx_t *rowptr = (idx_t *)crs->rowptr;
    const idx_t rstart = rowptr[0];

    if (rstart) {
        for (ptrdiff_t r = 0; r <= crs->lrows; r++) {
            rowptr[r] -= rstart;
        }
    }

    assert(crs->lnnz == rowptr[crs->lrows]);

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
    CATCH_MPI_ERROR(MPI_Type_size(crs->rowptr_type, &type_size));

    if (type_size != sizeof(idx_t)) {
        MPI_Abort(comm, -1);
    }

    CATCH_MPI_ERROR(MPI_Type_size(crs->colidx_type, &type_size));

    if (type_size != sizeof(idx_t)) {
        MPI_Abort(comm, -1);
    }

    int ret = ParMETIS_V3_PartKway(vtxdist,               // 0
                                   rowptr,                // 1
                                   (idx_t *)crs->colidx,  // 2
                                   vwgt,                  // 3
                                   adjwgt,                // 4
                                   &wgtflag,              // 5
                                   &numflag,              // 6
                                   &ncon,                 // 7
                                   &nparts,               // 8
                                   tpwgts,                // 9
                                   ubvec,                 // 10
                                   options,               // 11
                                   &edgecut,              // 12
                                   parts,                 // 13
                                   &comm);                // 14

    // Undo crs modification
    if (rstart) {
        for (ptrdiff_t r = 0; r <= crs->lrows; r++) {
            rowptr[r] += rstart;
        }
    }

    return ret;
}
