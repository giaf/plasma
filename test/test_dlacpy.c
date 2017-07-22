/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlacpy.c, normal z -> d, Sat Jul 22 09:25:13 2017
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

/***************************************************************************//**
 *
 * @brief Tests DLACPY
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_dlacpy(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_UPLO   ].used = true;
    param[PARAM_TRANSA ].used = true;
    param[PARAM_DIM    ].used = PARAM_USE_M | PARAM_USE_N;
    param[PARAM_PADA   ].used = true;
    param[PARAM_PADB   ].used = true;
    param[PARAM_NB     ].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters
    //================================================================
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);
    plasma_enum_t transa = plasma_trans_const(param[PARAM_TRANSA].c);

    int m = param[PARAM_DIM].dim.m;
    int n = param[PARAM_DIM].dim.n;

    int lda = imax(1, m + param[PARAM_PADA].i);
    int ldb = imax(1, m + param[PARAM_PADB].i);

    int    test = param[PARAM_TEST].c == 'y';
    double eps  = LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters
    //================================================================
    plasma_set(PlasmaTuning, PlasmaDisabled);
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays
    //================================================================
    double *A =
        (double*)malloc((size_t)lda*n*sizeof(double));
    assert(A != NULL);

    double *B =
        (double*)malloc((size_t)ldb*n*sizeof(double));
    assert(B != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_dlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    lapack_int mtrxLayout = LAPACK_COL_MAJOR;

    // LAPACKE_[zcds]lantr_work has a bug (returns 0)
    // in MKL <= 11.3.3 (at least). Fixed in LAPACK 3.6.1.
    // For now, zero out the opposite triangle and use lange.
    // @sa test_dtrmm

    // Enforce zeroes in general rectangle or upper or lower triangle of A
    switch (uplo) {
        case PlasmaLower:
            LAPACKE_dlaset_work(
                mtrxLayout, 'U', m-1, n-1, 0.0, 0.0, &A[m], lda);
            break;
        case PlasmaUpper:
            LAPACKE_dlaset_work(
                mtrxLayout, 'L', m-1, n-1, 0.0, 0.0, &A[1], lda);
            break;
    }
    // Zero out B
    int Bm = (transa == PlasmaGeneral ? m : n);
    int Bn = (transa == PlasmaGeneral ? n : m);
    LAPACKE_dlaset_work(
        mtrxLayout, 'G', Bm, Bn, 0.0, 0.0, B, ldb);

    //================================================================
    // Run and time PLASMA
    //================================================================
    plasma_time_t start = omp_get_wtime();

    retval = plasma_dlacpy(uplo, transa, m, n, A, lda, B, ldb);

    plasma_time_t stop = omp_get_wtime();

    param[PARAM_GFLOPS].d = 0.0;

    if (retval != PlasmaSuccess) {
        plasma_error("plasma_dlacpy() failed");
        param[PARAM_TIME].d    = 0.0;
        param[PARAM_ERROR].d   = 1.0;
        param[PARAM_SUCCESS].i = false;
        return;
    }
    else {
        param[PARAM_TIME].d = stop-start;
    }

    //================================================================
    // Test results by comparing to result of core_dlacpy function
    //================================================================
    if (test) {
        // Calculate relative error |op(A) - B|_F / |A|_F < 3*eps
        // Using 3*eps covers complex arithmetic

        // Calculate difference op(A)-B
        if (transa == PlasmaTrans) {
            for (int i=0; i<m; i++)
                for (int j=0; j<n; j++)
                    B[j + i*ldb] -= A[i + j*lda];
        } else if (transa == PlasmaConjTrans) {
            for (int i=0; i<m; i++)
                for (int j=0; j<n; j++)
                    B[j + i*ldb] -= (A[i + j*lda]);
        } else {
            for (int i=0; i<m; i++)
                for (int j=0; j<n; j++)
                    B[i + j*ldb] -= A[i + j*lda];
        }
        if (retval != PlasmaSuccess) {
            coreblas_error("LAPACKE_dlacpy_work() failed");
            param[PARAM_ERROR].d   = 1.0;
            param[PARAM_SUCCESS].i = false;
            return;
        }

        double work[1];

        // Calculate Frobenius norm of A
        double Anorm = LAPACKE_dlange_work(
                           mtrxLayout, 'F', m, n, A, lda, work);


        // Calculate Frobenius norm of op(A)-B
        double BnormDiff = LAPACKE_dlange_work(
                               mtrxLayout, 'F', Bm, Bn, B, ldb, work);

        // Calculate relative error |op(A)-B|_F / |A|_F
        double error = BnormDiff/Anorm;

        param[PARAM_ERROR].d   = error;
        param[PARAM_SUCCESS].i = error < 3*eps;
    }

    //================================================================
    // Free arrays
    //================================================================
    free(A);
    free(B);
}
