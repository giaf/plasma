/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zsyrk.c, normal z -> d, Thu Aug  8 17:24:59 2019
 *
 **/

#include <plasma_core_blas.h>
#include "plasma_types.h"
#include "core_lapack.h"
#include "blasfeo_d_aux.h"

/***************************************************************************//**
 *
 * @ingroup core_syrk
 *
 *  Performs one of the symmetric rank k operations
 *
 *    \f[ C = \alpha A \times A^T + \beta C, \f]
 *    or
 *    \f[ C = \alpha A^T \times A + \beta C, \f]
 *
 *  where alpha and beta are scalars, C is an n-by-n symmetric
 *  matrix, and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaUpper: Upper triangle of C is stored;
 *          - PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          - PlasmaNoTrans: \f[ C = \alpha A \times A^T + \beta C; \f]
 *          - PlasmaTrans:   \f[ C = \alpha A^T \times A + \beta C. \f]
 *
 * @param[in] n
 *          The order of the matrix C. n >= 0.
 *
 * @param[in] k
 *          If trans = PlasmaNoTrans, number of columns of the A matrix;
 *          if trans = PlasmaTrans, number of rows of the A matrix.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          A is an lda-by-ka matrix.
 *          If trans = PlasmaNoTrans, ka = k;
 *          if trans = PlasmaTrans,   ka = n.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          If trans = PlasmaNoTrans, lda >= max(1, n);
 *          if trans = PlasmaTrans,   lda >= max(1, k).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          C is an ldc-by-n matrix.
 *          On exit, the uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1, n).
 *
 ******************************************************************************/
__attribute__((weak))
void plasma_core_dsyrk_blasfeo(plasma_enum_t uplo, plasma_enum_t trans,
                int n, int k,
                double alpha, struct blasfeo_dmat *sA, int ai, int aj,
                double beta,  struct blasfeo_dmat *sC, int ci, int cj)
{
    // cblas_dsyrk(CblasColMajor,
    //             (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
    //             n, k,
    //             (alpha), A, lda,
    //             (beta),  C, ldc);
    blasfeo_dsyrk_ln(n, k, alpha, sA, ai, aj, sA, ai, aj, beta, sC, ci, cj, sC, ci, cj);
}

/******************************************************************************/
void plasma_core_omp_dsyrk_blasfeo(
    plasma_enum_t uplo, plasma_enum_t trans,
    int n, int k,
    double alpha, struct blasfeo_dmat *sA, int ai, int aj,
    double beta,  struct blasfeo_dmat *sC, int ci, int cj,
    plasma_sequence_t *sequence, plasma_request_t *request)
{
    int ak;
    if (trans == PlasmaNoTrans)
        ak = k;
    else
        ak = n;

	struct blasfeo_dmat sA2, sC2;
    sA2 = *sA;
	sC2 = *sC;

    double *A = sA->pA;
	int sda = sA->cn;
	double *C = sC->pA;
	int sdc = sC->cn;

    // #pragma omp task depend(in:A[0:lda*ak]) \
    //                  depend(inout:C[0:ldc*n])
    #pragma omp task depend(in:A[0:(sA->pm)*(sA->cn)]) \
                     depend(inout:C[0:(sC->pm)*(sC->cn)])
    {
        if (sequence->status == PlasmaSuccess)
            plasma_core_dsyrk_blasfeo(uplo, trans,
                       n, k,
                       alpha, &sA2, ai, aj,
                       beta,  &sC2, ci, cj);
    }
}
