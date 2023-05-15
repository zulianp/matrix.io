#include "matrixio_checks.h"

namespace detail {

    template <typename T>
    static inline bool crs_rowptr_is_sane(const ptrdiff_t nrows, const T *const rowptr) {
        bool sane = true;

        for (ptrdiff_t i = 0; i <= nrows; i++) {
            if (rowptr[i] < 0) {
                sane = false;
            }
        }

        for (ptrdiff_t i = 0; i < nrows; i++) {
            if (rowptr[i] > rowptr[i + 1]) {
                sane = false;
            }
        }

        return sane;
    }

    template <typename T>
    static inline bool crs_colidx_is_sane(const ptrdiff_t nnz, const ptrdiff_t ncols, const T *const colidx) {
        bool sane = true;
        for (ptrdiff_t i = 0; i < nnz; i++) {
            if (colidx[i] >= ncols) {
                sane = false;
            }
        }
        return sane;
    }

    template <typename RowPtr, typename ColIdx>
    static inline bool crs_graph_is_sane(const ptrdiff_t nrows,
                                         const ptrdiff_t ncols,
                                         const RowPtr *const rowptr,
                                         const ColIdx *const colidx) {
        return crs_rowptr_is_sane(nrows, rowptr) && crs_colidx_is_sane(rowptr[nrows], ncols, colidx);
    }

    // template <typename RowPtr, typename ColIdx>
    // static inline bool crs_graph_is_symmetric(const ptrdiff_t nrows,
    //                                           const ptrdiff_t ncols,
    //                                           const RowPtr *const rowptr,
    //                                           const ColIdx *const colidx) {

    //     const RowPtr first_row = rowptr[0];
    //     for(ptrdiff_t i = 0; i < nrows; i++) {
    //         const RowPtr ibegin = rowptr[i] - first_row;
    //         const RowPtr iextent = rowptr[i+1] - rowptr[i];
    //         const ColIdx *lcidx = &colidx[ibegin];

    //         for(RowPtr ik = 0; ik < iextent; ik++) {
    //             const ptrdiff_t j = lcidx[ik];

    //             const RowPtr ibegin = rowptr[i] - first_row;
    //             const RowPtr iextent = rowptr[i+1] - rowptr[i];
    //             const ColIdx *lcidx = &colidx[ibegin];
    //         }

    //     }

    // }

}  // namespace detail

extern "C" {

int crs_check_graph(const ptrdiff_t nrows,
                    const ptrdiff_t ncols,
                    const ptrdiff_t nnz,
                    MPI_Datatype rowptr_type,
                    const void *const rowptr,
                    MPI_Datatype colidx_type,
                    const void *const colidx) {
    bool rowptr_sane = false;
    if (rowptr_type == MPI_INT) {
        rowptr_sane = detail::crs_rowptr_is_sane<int>(nrows, (int *)rowptr);
    } else if (rowptr_type == MPI_LONG) {
        rowptr_sane = detail::crs_rowptr_is_sane<long>(nrows, (long *)rowptr);
    } else if (rowptr_type == MPI_UNSIGNED) {
        rowptr_sane = detail::crs_rowptr_is_sane<unsigned>(nrows, (unsigned *)rowptr);
    }

    bool colidx_sane = false;
    if (colidx_type == MPI_INT) {
        colidx_sane = detail::crs_colidx_is_sane<int>(nnz, ncols, (int *)colidx);
    } else if (colidx_type == MPI_LONG) {
        colidx_sane = detail::crs_colidx_is_sane<long>(nnz, ncols, (long *)colidx);
    } else if (colidx_type == MPI_UNSIGNED) {
        colidx_sane = detail::crs_colidx_is_sane<unsigned>(nnz, ncols, (unsigned *)colidx);
    }

    return !rowptr_sane + (!colidx_sane) * 2;
}

// int crs_graph_is_symmetric(const ptrdiff_t nrows,
//                            const ptrdiff_t ncols,
//                            const ptrdiff_t nnz,
//                            MPI_Datatype rowptr_type,
//                            const void *const rowptr,
//                            MPI_Datatype colidx_type,
//                            const void *const colidx) {}
}
