/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zcgesv.c, mixed zc -> ds, Sat Jul 22 09:25:26 2017
 *
 **/

#include "core_blas.h"
#include "core_lapack.h"
#include "flops.h"
#include "plasma.h"
#include "test.h"

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define REAL

#define A(i_, j_) A[(i_) + (size_t)lda*(j_)]

/***************************************************************************//**
 *
 * @brief Tests DSPOSV
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_dsgesv(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_DIM    ].used = PARAM_USE_N;
    param[PARAM_NRHS   ].used = true;
    param[PARAM_PADA   ].used = true;
    param[PARAM_PADB   ].used = true;
    param[PARAM_NB     ].used = true;
    param[PARAM_IB     ].used = true;
    param[PARAM_MTPF   ].used = true;
    param[PARAM_ITERSV ].used = true;
    param[PARAM_ZEROCOL].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters
    //================================================================
    int n    = param[PARAM_DIM].dim.n;
    int nrhs = param[PARAM_NRHS].i;
    int lda  = imax(1, n + param[PARAM_PADA].i);
    int ldb  = imax(1, n + param[PARAM_PADB].i);
    int ldx  = ldb;
    int ITER;

    int    test = param[PARAM_TEST].c == 'y';
    double tol  = param[PARAM_TOL].d * LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters
    //================================================================
    plasma_set(PlasmaTuning, PlasmaDisabled);
    plasma_set(PlasmaNb, param[PARAM_NB].i);
    plasma_set(PlasmaIb, param[PARAM_IB].i);
    plasma_set(PlasmaNumPanelThreads, param[PARAM_MTPF].i);

    //================================================================
    // Allocate and initialize arrays
    //================================================================
    double *A = (double *)malloc(
        (size_t)lda*n*sizeof(double));
    assert(A != NULL);

    int *ipiv = (int*)malloc((size_t)n*sizeof(int));
    assert(ipiv != NULL);

    double *B = (double *)malloc(
        (size_t)ldb*nrhs*sizeof(double));
    assert(B != NULL);

    double *X = (double *)malloc(
        (size_t)ldx*nrhs*sizeof(double));
    assert(X != NULL);

    // Initialize random A
    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_dlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    int zerocol = param[PARAM_ZEROCOL].i;
    if (zerocol >= 0 && zerocol < n)
        memset(&A[zerocol*lda], 0, n*sizeof(double));

    // Initialize B
    retval = LAPACKE_dlarnv(1, seed, (size_t)ldb*nrhs, B);
    assert(retval == 0);

    double *Aref = NULL;
    if (test) {
        Aref = (double *)malloc(
            (size_t)lda*n*sizeof(double));
        assert(Aref != NULL);
        memcpy(Aref, A, (size_t)lda*n*sizeof(double));
    }

    //================================================================
    // Run and time PLASMA
    //================================================================
    plasma_time_t start = omp_get_wtime();
    int plainfo = plasma_dsgesv(n, nrhs, A, lda, ipiv, B, ldb, X, ldx, &ITER);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;
    double flops = flops_dgetrf(n, n) + flops_dgetrs(n, nrhs);
    param[PARAM_ITERSV].i = ITER;
    param[PARAM_TIME].d   = time;
    param[PARAM_GFLOPS].d = flops / time / 1e9;

    //================================================================
    // Test results by checking the residual
    //
    //                      || B - AX ||_I
    //                --------------------------- < epsilon
    //                 || A ||_I * || X ||_I * N
    //
    //================================================================
    if (test) {
        if (plainfo == 0) {
            double alpha =  1.0;
            double beta  = -1.0;

            lapack_int mtrxLayout = LAPACK_COL_MAJOR;
            lapack_int mtrxNorm   = 'I';

            double *work = (double *)malloc(n*sizeof(double));
            assert(work != NULL);

            // Calculate infinite norms of matrices A_ref and X
            double Anorm = LAPACKE_dlange_work(mtrxLayout, mtrxNorm, n, n, Aref,
                                               lda, work);
            double Xnorm = LAPACKE_dlange_work(mtrxLayout, mtrxNorm, n, nrhs, X,
                                               ldx, work);

            // Calculate residual R = A*X-B, store result in B
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, nrhs, n,
                        (alpha), Aref, lda,
                                            X,    ldx,
                        (beta),  B,    ldb);

            // Calculate infinite norm of residual matrix R
            double Rnorm = LAPACKE_dlange_work(mtrxLayout, mtrxNorm, n, nrhs, B,
                                               ldb, work);
            // Calculate relative error
            double residual = Rnorm / ( n*Anorm*Xnorm );

            param[PARAM_ERROR].d   = residual;
            param[PARAM_SUCCESS].i = residual < tol;

            free(work);
        }
        else {
            int lapinfo = LAPACKE_dsgesv(
                              LAPACK_COL_MAJOR,
                              n, nrhs, A, lda, ipiv, B, ldb, X, ldx, &ITER);
            if (plainfo == lapinfo) {
                param[PARAM_ERROR].d = 0.0;
                param[PARAM_SUCCESS].i = 1;
            }
            else {
                param[PARAM_ERROR].d = INFINITY;
                param[PARAM_SUCCESS].i = 0;
            }
        }

        free(Aref);
    }

    //================================================================
    // Free arrays
    //================================================================
    free(A); free(ipiv); free(B); free(X);
}
