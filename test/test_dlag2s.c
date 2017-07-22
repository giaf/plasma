/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlag2c.c, mixed zc -> ds, Sat Jul 22 09:25:32 2017
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
 * @brief Tests DLAG2S
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_dlag2s(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_DIM    ].used = PARAM_USE_M | PARAM_USE_N;
    param[PARAM_PADA   ].used = true;
    param[PARAM_PADB   ].used = true;
    param[PARAM_NB     ].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters
    //================================================================
    int m = param[PARAM_DIM].dim.m;
    int n = param[PARAM_DIM].dim.n;

    int lda  = imax(1, m + param[PARAM_PADA].i);
    int ldas = imax(1, m + param[PARAM_PADB].i);

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

    float *As =
        (float*)malloc((size_t)ldas*n*sizeof(float));
    assert(As != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_dlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    float *AsRef = NULL;
    if (test) {
        AsRef = (float*)malloc(
            (size_t)ldas*n*sizeof(float));
        assert(AsRef != NULL);

        memcpy(AsRef, As, (size_t)ldas*n*sizeof(float));
    }

    //================================================================
    // Run and time PLASMA
    //================================================================
    plasma_time_t start = omp_get_wtime();

    retval = plasma_dlag2s(m, n, A, lda, As, ldas);

    plasma_time_t stop = omp_get_wtime();

    param[PARAM_GFLOPS].d  = 0.0;

    if (retval != PlasmaSuccess) {
        plasma_error("plasma_dlag2s() failed");
        param[PARAM_TIME].d    = 0.0;
        param[PARAM_ERROR].d   = 1.0;
        param[PARAM_SUCCESS].i = false;
        return;
    }
    else {
        param[PARAM_TIME].d = stop-start;
    }

    //================================================================
    // Test results by comparing to result of LAPACK function
    //================================================================
    if (test) {
        // Calculate relative error |As_ref - As|_F / |As_ref|_F < 3*eps
        // Using 3*eps covers complex arithmetic

        lapack_int mtrxLayout = LAPACK_COL_MAJOR;

        retval = LAPACKE_dlag2s_work(mtrxLayout, m, n, A, lda, AsRef, ldas);

        if (retval != PlasmaSuccess) {
            coreblas_error("LAPACKE_dlag2s_work() failed");
            param[PARAM_ERROR].d   = 1.0;
            param[PARAM_SUCCESS].i = false;
            return;
        }

        float work[1];

        // Calculate Frobenius norm of reference result As_ref
        double AsNormRef = LAPACKE_slange_work(mtrxLayout, 'F',
                                               ldas, n, AsRef, ldas, work);

        // Calculate difference As_ref-As
        float cmone = -1.0;
        cblas_saxpy((size_t)ldas*n, (cmone), As, 1, AsRef, 1);

        // Calculate Frobenius norm of As_ref-As
        double AsNormDiff = LAPACKE_slange_work(mtrxLayout, 'F',
                                                ldas, n, AsRef, ldas, work);

        // Calculate relative error |As_ref-As|_F / |As_ref|_F
        double error = AsNormDiff/AsNormRef;

        param[PARAM_ERROR].d   = error;
        param[PARAM_SUCCESS].i = error < 3*eps;
    }

    //================================================================
    // Free arrays
    //================================================================
    free(A);
    free(As);

    if (test)
        free(AsRef);
}
