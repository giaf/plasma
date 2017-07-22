/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlantr.c, normal z -> d, Sat Jul 22 09:25:31 2017
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
#include <math.h>

#include <omp.h>

#define REAL

/***************************************************************************//**
 *
 * @brief Tests DLANTR.
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_dlantr(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_NORM   ].used = true;
    param[PARAM_UPLO   ].used = true;
    param[PARAM_DIAG   ].used = true;
    param[PARAM_DIM    ].used = PARAM_USE_M | PARAM_USE_N;
    param[PARAM_PADA   ].used = true;
    param[PARAM_NB     ].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t norm = plasma_norm_const(param[PARAM_NORM].c);
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);
    plasma_enum_t diag = plasma_diag_const(param[PARAM_DIAG].c);

    int m = param[PARAM_DIM].dim.m;
    int n = param[PARAM_DIM].dim.n;

    int lda = imax(1, m + param[PARAM_PADA].i);

    int test = param[PARAM_TEST].c == 'y';
    double eps = LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaTuning, PlasmaDisabled);
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    double *A =
        (double*)malloc((size_t)lda*n*sizeof(double));
    assert(A != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_dlarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    double *Aref = NULL;
    double *work = NULL;
    if (test) {
        Aref = (double*)malloc(
            (size_t)lda*n*sizeof(double));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(double));

        work = (double*) malloc(imax(m, n)*sizeof(double));
        assert(work != NULL);
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    double value = plasma_dlantr(norm, uplo, diag, m, n, A, lda);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_dlange(m, n, norm) / time / 1e9;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        // LAPACKE_[ds]lantr_work has a bug (returns 0)
        // in MKL <= 11.3.3 (at least). Fixed in LAPACK 3.6.1.
        // For now, call LAPACK directly.
        // LAPACK_dlantr is a macro for correct name mangling (e.g.
        // adding _ at the end) of the Fortran symbol.
        // The macro is either defined in lapacke.h, or in the file
        // core_lapack_d.h for the use with MKL.
        char normc = lapack_const(norm);
        char uploc = lapack_const(uplo);
        char diagc = lapack_const(diag);
        double valueRef = LAPACK_dlantr(&normc, &uploc, &diagc,
                                        &m, &n, Aref, &lda, work);
        // double valueRef =
        //     LAPACKE_dlantr(LAPACK_COL_MAJOR,
        //                    lapack_const(norm), lapack_const(uplo),
        //                    lapack_const(diag),
        //                    m, n, Aref, lda);

        // Calculate relative error
        double error = fabs(value-valueRef);
        if (valueRef != 0)
            error /= valueRef;
        double tol = eps;
        double normalize = 1;
        switch (norm) {
            case PlasmaInfNorm:
                // Sum order on the line can differ
                normalize = n;
                break;

            case PlasmaOneNorm:
                // Sum order on the column can differ
                normalize = m;
                break;

            case PlasmaFrobeniusNorm:
                // Sum order on every element can differ
                normalize = m*n;
                break;
        }
        error /= normalize;
        param[PARAM_ERROR].d   = error;
        param[PARAM_SUCCESS].i = error < tol;
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    if (test) {
        free(Aref);
        free(work);
    }
}
