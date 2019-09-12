/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_ztrsm.c, normal z -> d, Thu Aug  8 17:24:56 2019
 *
 **/

#include <plasma_core_blas.h>
#include "plasma_types.h"
#include "core_lapack.h"
#include "blasfeo_d_aux.h"

/***************************************************************************//**
 *
 * @ingroup core_trsm
 *
 *  Solves one of the matrix equations
 *
 *    \f[ op( A )\times X  = \alpha B, \f] or
 *    \f[ X \times op( A ) = \alpha B, \f]
 *
 *  where op( A ) is one of:
 *    \f[ op( A ) = A,   \f]
 *    \f[ op( A ) = A^T, \f]
 *    \f[ op( A ) = A^T, \f]
 *
 *  alpha is a scalar, X and B are m-by-n matrices, and
 *  A is a unit or non-unit, upper or lower triangular matrix.
 *  The matrix X overwrites B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          - PlasmaLeft:  op(A)*X = B,
 *          - PlasmaRight: X*op(A) = B.
 *
 * @param[in] uplo
 *          - PlasmaUpper: A is upper triangular,
 *          - PlasmaLower: A is lower triangular.
 *
 * @param[in] transa
 *          - PlasmaNoTrans:   A is not transposed,
 *          - PlasmaTrans:     A is transposed,
 *          - PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          - PlasmaNonUnit: A has non-unit diagonal,
 *          - PlasmaUnit:    A has unit diagonal.
 *
 * @param[in] m
 *          The number of rows of the matrix B. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix B. n >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          The lda-by-ka triangular matrix,
 *          where ka = m if side = PlasmaLeft,
 *            and ka = n if side = PlasmaRight.
 *          If uplo = PlasmaUpper, the leading k-by-k upper triangular part
 *          of the array A contains the upper triangular matrix, and the
 *          strictly lower triangular part of A is not referenced.
 *          If uplo = PlasmaLower, the leading k-by-k lower triangular part
 *          of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced.
 *          If diag = PlasmaUnit, the diagonal elements of A are also not
 *          referenced and are assumed to be 1.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,k).
 *
 * @param[in,out] B
 *          On entry, the ldb-by-n right hand side matrix B.
 *          On exit, if return value = 0, the ldb-by-n solution matrix X.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void plasma_core_dtrsm_blasfeo(plasma_enum_t side, plasma_enum_t uplo,
                plasma_enum_t transa, plasma_enum_t diag,
                int m, int n,
                double alpha,  struct blasfeo_dmat *sA, int ai, int aj,
                               struct blasfeo_dmat *sB, int bi, int bj)
{
    // cblas_dtrsm(CblasColMajor,
    //             (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
    //             (CBLAS_TRANSPOSE)transa, (CBLAS_DIAG)diag,
    //             m, n,
    //             (alpha), A, lda,
    //                                 B, ldb);
	// TODO add checks for all unsupported variants !!!!!!!!!!!!!
    blasfeo_dtrsm_rltn(m, n, alpha, sA, ai, aj, sB, bi, bj, sB, bi, bj);
}

/******************************************************************************/
void plasma_core_omp_dtrsm_blasfeo(
    plasma_enum_t side, plasma_enum_t uplo,
    plasma_enum_t transa, plasma_enum_t diag,
    int m, int n,
    double alpha, const struct blasfeo_dmat *sA, int ai, int aj,
                                    struct blasfeo_dmat *sB, int bi, int bj,
    plasma_sequence_t *sequence, plasma_request_t *request)
{
    int ak;
    if (side == PlasmaLeft)
        ak = m;
    else
        ak = n;

    struct blasfeo_dmat sA2, sB2;
    sA2 = *sA;
	sB2 = *sB;

    double *A = sA->pA;
	int sda = sA->cn;
	double *B = sB->pA;
	int sdb = sB->cn;

    // #pragma omp task depend(in:A[0:lda*ak]) \
    //                  depend(inout:B[0:ldb*n])
    #pragma omp task depend(in:A[0:(sA->pm)*(sA->cn)]) \
                     depend(inout:B[0:(sB->pm)*(sB->cn)])
    {
        if (sequence->status == PlasmaSuccess)
            plasma_core_dtrsm_blasfeo(side, uplo,
                       transa, diag,
                       m, n,
                       alpha, &sA2, ai, aj,
                              &sB2, bi, bj);
    }
}
