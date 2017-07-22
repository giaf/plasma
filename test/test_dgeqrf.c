/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zgeqrf.c, normal z -> d, Sat Jul 22 09:25:40 2017
 *
 **/
#include "test.h"
#include "flops.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define REAL

/***************************************************************************//**
 *
 * @brief Tests DGEQRF.
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_dgeqrf(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_DIM    ].used = PARAM_USE_M | PARAM_USE_N;
    param[PARAM_PADA   ].used = true;
    param[PARAM_NB     ].used = true;
    param[PARAM_IB     ].used = true;
    param[PARAM_HMODE  ].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters.
    //================================================================
    int m = param[PARAM_DIM].dim.m;
    int n = param[PARAM_DIM].dim.n;

    int lda = imax(1, m + param[PARAM_PADA].i);

    int test = param[PARAM_TEST].c == 'y';
    double tol = param[PARAM_TOL].d * LAPACKE_dlamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaTuning, PlasmaDisabled);
    plasma_set(PlasmaNb, param[PARAM_NB].i);
    plasma_set(PlasmaIb, param[PARAM_IB].i);
    if (param[PARAM_HMODE].c == 't')
        plasma_set(PlasmaHouseholderMode, PlasmaTreeHouseholder);
    else
        plasma_set(PlasmaHouseholderMode, PlasmaFlatHouseholder);

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
    if (test) {
        Aref = (double*)malloc(
            (size_t)lda*n*sizeof(double));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(double));
    }

    //================================================================
    // Prepare the descriptor for matrix T.
    //================================================================
    plasma_desc_t T;

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    plasma_dgeqrf(m, n, A, lda, &T);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_dgeqrf(m, n) / time / 1e9;

    //=================================================================
    // Test results by checking orthogonality of Q and precision of Q*R
    //=================================================================
    if (test) {
        // Check the orthogonality of Q
        int minmn = imin(m, n);

        // Allocate space for Q.
        int ldq = m;
        double *Q =
            (double *)malloc((size_t)ldq*minmn*
                                         sizeof(double));

        // Build Q.
        plasma_dorgqr(m, minmn, minmn, A, lda, T, Q, ldq);

        // Build the identity matrix
        double *Id =
            (double *) malloc((size_t)minmn*minmn*
                                          sizeof(double));
        LAPACKE_dlaset_work(LAPACK_COL_MAJOR, 'g', minmn, minmn,
                            0.0, 1.0, Id, minmn);

        // Perform Id - Q^T * Q
        cblas_dsyrk(CblasColMajor, CblasUpper, CblasConjTrans, minmn, m,
                    -1.0, Q, ldq, 1.0, Id, minmn);

        // work array of size m is needed for computing L_oo norm
        double *work = (double *) malloc((size_t)m*sizeof(double));

        // |Id - Q^T * Q|_oo
        double ortho = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, 'I', 'u',
                                           minmn, Id, minmn, work);

        // normalize the result
        // |Id - Q^T * Q|_oo / n
        ortho /= minmn;

        free(Q);
        free(Id);

        // Check the accuracy of A - Q * R
        // LAPACK version does not construct Q, it uses only application of it

        // Extract the R.
        double *R =
            (double *)malloc((size_t)m*n*
                                         sizeof(double));
        LAPACKE_dlaset_work(LAPACK_COL_MAJOR, 'l', m, n,
                            0.0, 0.0, R, m);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'u', m, n, A, lda, R, m);

        // Compute Q * R.
        plasma_dormqr(PlasmaLeft, PlasmaNoTrans, m, n, minmn, A, lda, T,
                      R, m);

        // Compute the difference.
        // R = A - Q*R
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                R[j*m+i] = Aref[j*lda+i] - R[j*m+i];

        // |A|_oo
        double normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', m, n,
                                           Aref, lda, work);

        // |A - Q*R|_oo
        double error = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', m, n,
                                               R, m, work);

        // normalize the result
        // |A-QR|_oo / (|A|_oo * n)
        error /= (normA * n);

        param[PARAM_ERROR].d = error;
        param[PARAM_ORTHO].d = ortho;
        param[PARAM_SUCCESS].i = (error < tol && ortho < tol);

        free(work);
        free(R);
    }

    //================================================================
    // Free arrays.
    //================================================================
    plasma_desc_destroy(&T);
    free(A);
    if (test)
        free(Aref);
}
