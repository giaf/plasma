/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgetri_aux.c, normal z -> c, Sat Jul 22 09:25:44 2017
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

#define COMPLEX

#define A(i_, j_) A[(i_) + (size_t)lda*(j_)]

/***************************************************************************//**
 *
 * @brief Tests CTRSM.
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_cgetri_aux(param_value_t param[], bool run)
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
    plasma_complex32_t *A =
        (plasma_complex32_t*)malloc((size_t)lda*n*sizeof(plasma_complex32_t));
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
    retval = LAPACKE_clarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);
    LAPACKE_cgetrf(CblasColMajor, n, n, A, lda, ipiv);
    plasma_complex32_t *L = NULL;
    plasma_complex32_t *U = NULL;
    if (test) {
        L = (plasma_complex32_t*)malloc((size_t)n*n*sizeof(plasma_complex32_t));
        assert(L != NULL);
        U = (plasma_complex32_t*)malloc((size_t)n*n*sizeof(plasma_complex32_t));
        assert(U != NULL);

        LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'F', n, n, A, lda, L, n);
        LAPACKE_clacpy_work(LAPACK_COL_MAJOR, 'F', n, n, A, lda, U, n);

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

    plasma_cgetri_aux( n, A, lda );

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_ctrsm(PlasmaRight, n, n) / time / 1e9;

    //================================================================
    // Test results by checking the residual
    // ||A*L - B|| / (||A||*||L||)
    //================================================================
    if (test) {
        plasma_complex32_t zone  =  1.0;
        plasma_complex32_t zmone = -1.0;
        float work[1];
        float Anorm = LAPACKE_clange_work(
                           LAPACK_COL_MAJOR, 'F', n, n, A, lda, work);
        float Lnorm = LAPACKE_clange_work(
                           LAPACK_COL_MAJOR, 'F', n, n, L, n, work);

        // A*L - U
        cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
                    CBLAS_SADDR(zone),  A, lda,
                                        L, n,
                    CBLAS_SADDR(zmone), U, n);

        float error = LAPACKE_clange_work(
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
