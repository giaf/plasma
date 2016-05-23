/**
 *
 * @file core_zsyrk.c
 *
 *  PLASMA core_blas kernel.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley, Univ. of Colorado Denver and
 *  Univ. of Manchester.
 *
 * @version 3.0.0
 * @author Pedro V. Lara
 * @date 2016-05-11
 * @precisions normal z -> c d s
 *
 **/

#include "core_blas.h"
#include "plasma_types.h"

#ifdef PLASMA_WITH_MKL
    #include <mkl.h>
#else
    #include <cblas.h>
    #include <lapacke.h>
#endif

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  Performs one of the symmetric rank k operations
 *
 *    \f[ C = \alpha [ op( A ) \times conjg( op( A )' )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = conjg( X' )
 *
 *  where alpha and beta are real scalars, C is an n-by-n symmetric
 *  matrix and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          - PlasmaUpper: Upper triangle of C is stored;
 *          - PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or conjugate transposed:
 *          - PlasmaNoTrans:   A is not transposed;
 *          - PlasmaTrans  :   A is transposed.
 *
 * @param[in] n
 *          The number of n specifies the order of the matrix C. n must be at least zero.
 *
 * @param[in] k
 *          The number of k specifies the number of columns of the matrix op( A ).
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a lda-by-ka matrix, where ka is k when trans = PlasmaNoTrans,
 *          and is n otherwise.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda must be at least
 *          max( 1, n ) if trans == PlasmaNoTrans, otherwise lda must
 *          be at least max( 1, k ).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a ldc-by-n matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max( 1, n ).
 *
 *******************************************************************************/
void CORE_zsyrk(
    PLASMA_enum uplo, PLASMA_enum trans, 
    int n, int k,
    PLASMA_Complex64_t alpha, 
    const PLASMA_Complex64_t *A, int lda,
    PLASMA_Complex64_t beta,  
    PLASMA_Complex64_t *C, int ldc)
{

    cblas_zsyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        n, k,
        CBLAS_SADDR(alpha), 
	    A, lda,
        CBLAS_SADDR(beta), 
	    C, ldc);

}

/******************************************************************************/

void CORE_OMP_zsyrk(
    PLASMA_enum uplo, PLASMA_enum trans, 
    int n, int k,
    PLASMA_Complex64_t alpha, 
    const PLASMA_Complex64_t *A, int lda,
    PLASMA_Complex64_t beta,  
    PLASMA_Complex64_t *C, int ldc)
{
    // omp depends assume lda == n or k, and ldc == n,
    // depending on transposes
#pragma omp task depend(in:A[0:n*k]) depend(inout:C[0:n*n])
    CORE_zsyrk(
        uplo, trans,
        n, k,
        alpha, A, lda,
        beta, C, ldc);

}
