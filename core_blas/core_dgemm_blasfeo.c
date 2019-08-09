/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from core_blas/core_zgemm.c, normal z -> d, Thu Aug  8 10:20:04 2019
 *
 **/

#include <plasma_core_blas.h>
#include "plasma_types.h"
#include "core_lapack.h"

#include "blasfeo_d_aux.h"
//#include "blasfeo_d_blas.h"

/***************************************************************************//**
 *
 * @ingroup core_gemm
 *
 *  Performs one of the matrix-matrix operations
 *
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C, \f]
 *
 *  where op( X ) is one of:
 *    \f[ op( X ) = X,   \f]
 *    \f[ op( X ) = X^T, \f]
 *    \f[ op( X ) = X^T, \f]
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m-by-k matrix, op( B ) a k-by-n matrix and C an m-by-n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transa
 *          - PlasmaNoTrans:   A is not transposed,
 *          - PlasmaTrans:     A is transposed,
 *          - PlasmaConjTrans: A is conjugate transposed.
 *
 * @param[in] transb
 *          - PlasmaNoTrans:   B is not transposed,
 *          - PlasmaTrans:     B is transposed,
 *          - PlasmaConjTrans: B is conjugate transposed.
 *
 * @param[in] m
 *          The number of rows of the matrix op( A ) and of the matrix C.
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix op( B ) and of the matrix C.
 *          n >= 0.
 *
 * @param[in] k
 *          The number of columns of the matrix op( A ) and the number of rows
 *          of the matrix op( B ). k >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 *          An lda-by-ka matrix, where ka is k when transa = PlasmaNoTrans,
 *          and is m otherwise.
 *
 * @param[in] lda
 *          The leading dimension of the array A.
 *          When transa = PlasmaNoTrans, lda >= max(1,m),
 *          otherwise, lda >= max(1,k).
 *
 * @param[in] B
 *          An ldb-by-kb matrix, where kb is n when transb = PlasmaNoTrans,
 *          and is k otherwise.
 *
 * @param[in] ldb
 *          The leading dimension of the array B.
 *          When transb = PlasmaNoTrans, ldb >= max(1,k),
 *          otherwise, ldb >= max(1,n).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 *          An ldc-by-n matrix. On exit, the array is overwritten by the m-by-n
 *          matrix ( alpha*op( A )*op( B ) + beta*C ).
 *
 * @param[in] ldc
 *          The leading dimension of the array C. ldc >= max(1,m).
 *
 ******************************************************************************/
__attribute__((weak))
void plasma_core_dgemm_blasfeo(plasma_enum_t transa, plasma_enum_t transb,
                int m, int n, int k,
//                double alpha, const double *A, int lda,
//                                          const double *B, int ldb,
//                double beta,        double *C, int ldc)
                double alpha, const struct blasfeo_dmat *sA, int ai, int aj,
                                          const struct blasfeo_dmat *sB, int bi, int bj,
                double beta,        struct blasfeo_dmat *sC, int ci, int cj)
{
#if HAVE_BLASFEO_API
	// TODO assume double precision !!!
	// TODO assume nn version !!!
//printf("\n%d %d %d %d %d %d\n", m, n, k, lda, ldb, ldc);
//printf("\nhere\n");
//	struct blasfeo_dmat sA, sB, sC;
//	blasfeo_create_dmat(m, k, &sA, A);
//	sA.cn = lda;
//	blasfeo_create_dmat(k, n, &sB, B);
//	sB.cn = ldb;
//	blasfeo_create_dmat(m, n, &sC, C);
//	sC.cn = ldc;
//printf("\nafter create\n");
//	blasfeo_print_dmat(m, k, &sA, 0, 0);
//	blasfeo_print_dmat(k, n, &sB, 0, 0);
//	blasfeo_print_dmat(m, n, &sC, 0, 0);
//printf("\nafter print\n");
//	blasfeo_dgemm_nn(m, n, k, alpha, &sA, 0, 0, &sB, 0, 0, beta, &sC, 0, 0, &sC, 0, 0);
	// TODO fix 'n'
	if(transa=='n' && transb=='n')
		{
		blasfeo_dgemm_nn(m, n, k, alpha, sA, ai, aj, sB, bi, bj, beta, sC, ci, cj, sC, ci, cj);
		}
	else
		{
		printf("\nimplement the others too\n");
		}
//printf("\nafter dgemm\n");
//	blasfeo_print_dmat(m, n, &sC, 0, 0);
#else
    cblas_dgemm(CblasColMajor,
                (CBLAS_TRANSPOSE)transa, (CBLAS_TRANSPOSE)transb,
                m, n, k,
                (alpha), A, lda,
                                    B, ldb,
                (beta),  C, ldc);
#endif
}

/******************************************************************************/
void plasma_core_omp_dgemm_blasfeo(
    plasma_enum_t transa, plasma_enum_t transb,
    int m, int n, int k,
//    double alpha, const double *A, int lda,
    double alpha, const struct blasfeo_dmat *sA, int ai, int aj,
//                              const double *B, int ldb,
                              const struct blasfeo_dmat *sB, int bi, int bj,
//    double beta,        double *C, int ldc,
    double beta,        struct blasfeo_dmat *sC, int ci, int cj,
    plasma_sequence_t *sequence, plasma_request_t *request)
{
	struct blasfeo_dmat sA2, sB2, sC2;
//	blasfeo_create_dmat(m, k, &sA, A);
//	blasfeo_create_dmat(k, n, &sB, B);
//	blasfeo_create_dmat(m, n, &sC2, C);
//	sC2.cn = sdc;
	// make a local copy of the structure, such that the orignal one can be safely destroyed when out of scope
	sA2 = *sA;
	sB2 = *sB;
	sC2 = *sC;

    int ak, am;
    if (transa == PlasmaNoTrans)
		{
        am = m;
//        ak = k;
		}
    else
		{
        am = k;
//        ak = m;
		}
	am = (am+D_PS-1)/D_PS*D_PS;

    int bk;
    if (transb == PlasmaNoTrans)
		{
        bk = k;
//        bk = n;
		}
    else
		{
        bk = n;
//        bk = k;
		}
	bk = (bk+D_PS-1)/D_PS*D_PS;

	double *A = sA->pA;
	int sda = sA->cn;
	double *B = sB->pA;
	int sdb = sB->cn;
	double *C = sC->pA;
	int sdc = sC->cn;

//    #pragma omp task depend(in:A[0:lda*ak]) \
//                     depend(in:B[0:ldb*bk]) \
//                     depend(inout:C[0:ldc*n])
	// TODO this should be checked in case of ai,aj!=0
    #pragma omp task depend(in:A[0:am*sda]) \
                     depend(in:B[0:bk*sdb]) \
                     depend(inout:C[0:m*sdc])
    {
        if (sequence->status == PlasmaSuccess)
            plasma_core_dgemm_blasfeo(transa, transb,
                       m, n, k,
//                       alpha, A, lda,
//                              B, ldb,
//                       beta,  C, ldc);
                       alpha, &sA2, ai, aj,
                              &sB2, bi, bj,
                       beta,  &sC2, ci, cj);
    }
}
