#include "matrixio_crs.h"

#include "utils.h"

#include "parmetis.h"
#include <mpi.h>

#include <stdlib.h>
#include <stddef.h>


#define MPI_IDX_T MPI_INT

int decompose(MPI_Comm comm, crs_t*crs, int*const parts) {
	int rank, size;
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &size);

	idx_t *vtxdist = (idx_t *)malloc((size + 1) * sizeof(idx_t));
	// memset(vtxdist, 0, (size + 1) * sizeof(idx_t));

	idx_t rowoffset = crs->rowoffset;
	MPI_Allgather(
	    &rowoffset, 1, MPI_IDX_T, vtxdist, 1, MPI_IDX_T, comm);


	idx_t ncon = 1;
	idx_t *vwgt = 0;
	// idx_t* vsize = 0;
	idx_t *adjwgt = 0;
	idx_t wgtflag = 0;
	idx_t numflag = 0;
	idx_t nparts = size;

	real_t *tpwgts = (real_t*) malloc(size * sizeof(real_t));
	for(int r = 0; r < size; r++) {
		tpwgts[r] = 1./size;
	}


	real_t ubvec[1] = {1.05};
	idx_t options[3] = {0};
	// idx_t objval = -1;
	idx_t edgecut;

	int type_size;
	CATCH_MPI_ERROR(MPI_Type_size(crs->rowptr_type, &type_size));

	if(type_size != sizeof(idx_t)) {
		MPI_Abort(comm, -1);
	}

	CATCH_MPI_ERROR(MPI_Type_size(crs->colidx_type, &type_size));
	
	if(type_size != sizeof(idx_t)) {
		MPI_Abort(comm, -1);
	}

	int ret = ParMETIS_V3_PartKway(vtxdist,     // 0
	                               (idx_t*) crs->rowptr,  // 1
	                               (idx_t*)	crs->colidx,  // 2
	                               vwgt,        // 3
	                               adjwgt,      // 4
	                               &wgtflag,    // 5
	                               &numflag,    // 6
	                               &ncon,       // 7
	                               &nparts,     // 8
	                               tpwgts,  // 9
	                               ubvec,       // 10
	                               options,     // 11
	                               &edgecut,    // 12
	                               parts,       // 13
	                               &comm);      // 14

	return ret;
}
