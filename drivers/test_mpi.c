
#include <mpi.h>
#include <stdlib.h>

// mpicc test_mpi.c test_mpi
// ./test_mpi
// python3 -c "import numpy as np; print(np.fromfile('output_dump.raw', dtype=np.float64))"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    const char *path = "output_dump.raw";

    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    long nlocal = 10;
    long nglobal = 10 * size;
    double *data = malloc(nlocal * sizeof(double));

    for(long i = 0; i < nlocal; i++) {
    	data[i] = rank * size + i;
    }

    MPI_Datatype type = MPI_DOUBLE;

    MPI_Status status;
    MPI_Offset nbytes;
    MPI_File file;
    int type_size;

    MPI_Type_size(type, &type_size);
    nbytes = nglobal * type_size;

    MPI_File_open(comm, path, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &file);
    MPI_File_set_size(file, nbytes);

    long lnl = nlocal;
    long offset = 0;

    if (size > 1) {
        MPI_Exscan(&lnl, &offset, 1, MPI_LONG, MPI_SUM, comm);
    }

    MPI_File_write_at_all(file, offset * type_size, data, nlocal, type, &status);
    MPI_File_close(&file);

    free(data);

    return MPI_Finalize();
}
