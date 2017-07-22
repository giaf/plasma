/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgetri_aux.c, normal z -> s, Sat Jul 22 09:25:45 2017
 *
 **/
#include "test.h"
#include "flops.h"
#include "core_blas.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define REAL

#define A(i_, j_) A[(i_) + (size_t)lda*(j_)]

/***************************************************************************//**
 *
 * @brief Tests STRSM.
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_sgetri_aux(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_DIM    ].used = PARAM_USE_N;
    param[PARAM_PADA   ].used = true;
    param[PARAM_NB     ].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters.
    //================================================================
    int n = param[PARAM_DIM].dim.n;
    int lda = imax(1, n + param[PARAM_PADA].i);

    int test = param[PARAM_TEST].c == 'y';
    float tol = param[PARAM_TOL].d * LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaTuning, PlasmaDisabled);
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    float *A =
        (float*)malloc((size_t)lda*n*sizeof(float));
    assert(A != NULL);

    int *ipiv;
    ipiv = (int*)malloc((size_t)n*sizeof(int));
    assert(ipiv != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;

    //=================================================================
    // Initialize the matrices.
    // Factor A into LU to get well-conditioned triangular matrices.
    // Use L for unit triangle, and U for non-unit triangle,
    // transposing as necessary.
    // (There is some danger, as L^T or U^T may be much worse conditioned
    // than L or U, but in practice it seems okay.
    // See Higham, Accuracy and Stability of Numerical Algorithms, ch 8.)
    //=================================================================
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);
    LAPACKE_sgetrf(CblasColMajor, n, n, A, lda, ipiv);
    float *L = NULL;
    float *U = NULL;
    if (test) {
        L = (float*)malloc((size_t)n*n*sizeof(float));
        assert(L != NULL);
        U = (float*)malloc((size_t)n*n*sizeof(float));
        assert(U != NULL);

        LAPACKE_slacpy_work(LAPACK_COL_MAJOR, 'F', n, n, A, lda, L, n);
        LAPACKE_slacpy_work(LAPACK_COL_MAJOR, 'F', n, n, A, lda, U, n);

        for (int j = 0; j < n; j++) {
            L[j + j*n] = 1.0;
            for (int i = 0; i < j; i++) L[i + j*n] = 0.0;
            for (int i = j+1; i < n; i++) U[i + j*n] = 0.0;
        }
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();

    plasma_sgetri_aux( n, A, lda );

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_strsm(PlasmaRight, n, n) / time / 1e9;

    //================================================================
    // Test results by checking the residual
    // ||A*L - B|| / (||A||*||L||)
    //================================================================
    if (test) {
        float zone  =  1.0;
        float zmone = -1.0;
        float work[1];
        float Anorm = LAPACKE_slange_work(
                           LAPACK_COL_MAJOR, 'F', n, n, A, lda, work);
        float Lnorm = LAPACKE_slange_work(
                           LAPACK_COL_MAJOR, 'F', n, n, L, n, work);

        // A*L - U
        cblas_sgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
                    (zone),  A, lda,
                                        L, n,
                    (zmone), U, n);

        float error = LAPACKE_slange_work(
                           LAPACK_COL_MAJOR, 'F', n, n, U, n, work);
        if (Anorm*Lnorm != 0)
            error /= (Anorm * Lnorm);

        param[PARAM_ERROR].d = error;
        param[PARAM_SUCCESS].i = error < tol;
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(ipiv);
    if (test) {
        free(L); free(U);
    }
}
