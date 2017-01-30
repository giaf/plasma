/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> c d s
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "core_lapack.h"

#define COMPLEX
#define A(m,n) ((plasma_complex64_t*)plasma_tile_addr(A, m, n))
#define W2(j)  ((plasma_complex64_t*)plasma_tile_addr(W, j+A.mt, 0))   // 2mt x nb*nb


/***************************************************************************//**
 *
 * @ingroup core_laswp_sym
 *
 *  Applies symmetric pivoting.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaGeneral: entire A,
 *          - PlasmaUpper:   upper triangle,
 *          - PlasmaLower:   lower triangle.
 *
 * @param[in] A
 *          The matrix to be pivoted.
 *
 *
 *
 ******************************************************************************/
void core_zlaswp_sym(int uplo, plasma_desc_t A, int k1, int k2, const int *ipiv, int incx)
{
    if (uplo == PlasmaLower) {
        if (incx > 0) {
            for (int i = k1-1; i <= k2-1; i += incx) {
                if (ipiv[i]-1 != i) {
                    int p1 = i;
                    int p2 = ipiv[i]-1;

                    int i1 = p1%A.mb;
                    int i2 = p2%A.mb;
                    int m1 = p1/A.mb;
                    int m2 = p2/A.mb;
                    int lda1 = plasma_tile_mmain(A, m1);
                    int lda2 = plasma_tile_mmain(A, m2);


                    int i1p1 = (p1+1)%A.mb;
                    int i2p1 = (p2+1)%A.mb;
                    int m1p1 = (p1+1)/A.mb;
                    int m2p1 = (p2+1)/A.mb;
                    int lda1p1 = plasma_tile_mmain(A, m1p1);
                    int lda2p1 = plasma_tile_mmain(A, m2p1);

                    // swap rows of previous column (assuming (k1,k2) stay within a tile)
                    if (i > k1-1) {
                        cblas_zswap(i-(k1-1),
                                    A(m1, m1) + i1, lda1,
                                    A(m2, m1) + i2, lda2);
                    }

                    // swap columns p1 and p2
                    int mvam = plasma_tile_mview(A, m2p1);
                    if (mvam > i2+1) {
                        // between first tiles A(p2,p1) and A(p2,p2)
                        cblas_zswap(mvam-(i2+1),
                                    A(m2p1, m1) + i2p1 + i1*lda2p1, 1,
                                    A(m2p1, m2) + i2p1 + i2*lda2p1, 1);
                    }
                    for (int k = m2+1; k < A.mt; k++) {
                        int mvak = plasma_tile_mview(A, k);
                        int ldak = plasma_tile_mmain(A, k);
                        cblas_zswap(mvak,
                                    A(k, m1) + i1*ldak, 1,
                                    A(k, m2) + i2*ldak, 1);
                    }

                    // sym swap 
                    mvam = plasma_tile_mview(A, m1);
                    if (imin(mvam,p2-(k1-1)) > i1+1) {
                        #ifdef COMPLEX
                        LAPACKE_zlacgv_work(imin(mvam,p2-(k1-1))-(i1+1), A(m1p1, m1) + i1p1 + i1*lda1p1, 1);
                        LAPACKE_zlacgv_work(imin(mvam,p2-(k1-1))-(i1+1), A(m2, m1p1) + i2 + i1p1*lda2, lda2);
                        #endif
                        cblas_zswap(imin(mvam,p2-(k1-1))-(i1+1),
                                    A(m1p1, m1) + i1p1 + i1*lda1p1, 1,
                                    A(m2, m1p1) + i2 + i1p1*lda2, lda2);
                    }
                    for (int k = m1+1; k <= m2; k++) {
                        int mvak = plasma_tile_mview(A, k);
                        int ldak = plasma_tile_mmain(A, k);
                        #ifdef COMPLEX
                        LAPACKE_zlacgv_work(imin(mvak, (p2-1)-k*A.mb+1), A(k, m1) +  i1*ldak, 1);
                        LAPACKE_zlacgv_work(imin(mvak, (p2-1)-k*A.mb+1), A(m2, k) +  i2, lda2);
                        #endif
                        cblas_zswap(imin(mvak, (p2-1)-k*A.mb+1),
                                    A(k, m1) +  i1*ldak, 1,
                                    A(m2, k) +  i2, lda2);
                    }
                    #ifdef COMPLEX
                    LAPACKE_zlacgv_work(1, A(m2, m1) +  i2 + i1*lda2, 1);
                    #endif

                    // swap diagonal
                    cblas_zswap(1,
                                A(m1, m1) + i1 + i1*lda1, lda1,
                                A(m2, m2) + i2 + i2*lda2, lda2);

                }
            }
        }
    }
}

/******************************************************************************/
void core_omp_zlaswp_sym(int uplo,
                         plasma_desc_t A, int k1, int k2, const int *ipiv, int incx,
                         plasma_sequence_t *sequence, plasma_request_t *request)
{
    {
        if (sequence->status == PlasmaSuccess) {
            core_zlaswp_sym(uplo, A, k1, k2, ipiv, incx);
        }
    }
}