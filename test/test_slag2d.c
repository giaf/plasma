/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_clag2z.c, mixed zc -> ds, Sat Jul 22 09:25:41 2017
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
 * @brief Tests SLAG2D
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_slag2d(param_value_t param[], bool run)
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

    int ldas = imax(1, m + param[PARAM_PADA].i);
    int lda  = imax(1, m + param[PARAM_PADB].i);

    int    test = param[PARAM_TEST].c == 'y';
    double eps  = LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays
    //================================================================
    float *As =
        (float*)malloc((size_t)ldas*n*sizeof(float));
    assert(As != NULL);

    double *A =
        (double*)malloc((size_t)lda*n*sizeof(double));
    assert(A != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_slarnv(1, seed, (size_t)ldas*n, As);
    assert(retval == 0);

    double *Aref = NULL;
    if (test) {
        Aref = (double*)malloc(
            (size_t)ldas*n*sizeof(double));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(double));
    }

    //================================================================
    // Run and time PLASMA
    //================================================================
    plasma_time_t start = omp_get_wtime();

    retval = plasma_slag2d(m, n, As, ldas, A, lda);

    plasma_time_t stop = omp_get_wtime();

    param[PARAM_GFLOPS].d = 0.0;

    if (retval != PlasmaSuccess) {
        plasma_error("plasma_slag2d() failed");
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
        // Calculate relative error |A_ref - A|_F / |A_ref|_F < 3*eps
        // Using 3*eps covers complex arithmetic

        lapack_int mtrxLayout = LAPACK_COL_MAJOR;

        retval = LAPACKE_slag2d_work(mtrxLayout, m, n, As, ldas, Aref, lda);

        if (retval != PlasmaSuccess) {
            coreblas_error("LAPACKE_slag2d_work() failed");
            param[PARAM_ERROR].d   = 1.0;
            param[PARAM_SUCCESS].i = false;
            return;
        }

        double work[1];

        // Calculate Frobenius norm of reference result A_ref
        double AnormRef = LAPACKE_dlange_work(mtrxLayout, 'F',
                                              lda, n, Aref, lda, work);

        // Calculate difference A_ref-A
        double zmone = -1.0;
        cblas_daxpy((size_t)lda*n, (zmone), A, 1, Aref, 1);

        // Calculate Frobenius norm of A_ref-A
        double AnormDiff = LAPACKE_dlange_work(mtrxLayout, 'F',
                                               lda, n, Aref, lda, work);

        // Calculate relative error |A_ref-A|_F / |A_ref|_F
        double error = AnormDiff/AnormRef;

        param[PARAM_ERROR].d   = error;
        param[PARAM_SUCCESS].i = error < 3*eps;
    }

    //================================================================
    // Free arrays
    //================================================================
    free(As);
    free(A);

    if (test)
        free(Aref);
}
