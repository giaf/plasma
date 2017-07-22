/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zhetrf.c, normal z -> c, Sat Jul 22 09:25:25 2017
 *
 **/
#include "test.h"
#include "flops.h"
#include "core_blas.h"
#include "core_lapack.h"
#include "plasma.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <omp.h>

#define A(i_, j_) A[(i_) + (size_t)lda*(j_)]

/***************************************************************************//**
 *
 * @brief Tests CHETRF.
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_chetrf(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_UPLO   ].used = true;
    param[PARAM_DIM    ].used = PARAM_USE_N;
    param[PARAM_PADA   ].used = true;
    param[PARAM_NRHS   ].used = true;
    param[PARAM_NB     ].used = true;
    param[PARAM_IB     ].used = true;
    param[PARAM_MTPF   ].used = true;
    param[PARAM_ZEROCOL].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t uplo = plasma_uplo_const(param[PARAM_UPLO].c);

    int n = param[PARAM_DIM].dim.n;
    int nb = param[PARAM_NB].i;

    int lda = imax(1, n + param[PARAM_PADA].i);
    // band matrix A in skewed LAPACK storage
    int kut  = (nb+nb+nb-1)/nb; // # of tiles in upper band (not including diagonal)
    int klt  = (nb+nb-1)/nb;    // # of tiles in lower band (not including diagonal)
    int ldt  = (kut+klt+1)*nb;  // since we use cgetrf on panel, we pivot back within panel.

    int test = param[PARAM_TEST].c == 'y';
    float tol = param[PARAM_TOL].d * LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaTuning, PlasmaDisabled);
    plasma_set(PlasmaNb, param[PARAM_NB].i);
    plasma_set(PlasmaIb, param[PARAM_IB].i);
    plasma_set(PlasmaNumPanelThreads, param[PARAM_MTPF].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    plasma_complex32_t *A =
        (plasma_complex32_t*)calloc((size_t)lda*n, sizeof(plasma_complex32_t));
    assert(A != NULL);
    plasma_complex32_t *T =
        (plasma_complex32_t*)calloc((size_t)ldt*n, sizeof(plasma_complex32_t));
    assert(T != NULL);

    int *ipiv = (int*)calloc((size_t)n, sizeof(int));
    int *ipiv2 = (int*)calloc((size_t)n, sizeof(int));
    assert(ipiv != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_clarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);

    //================================================================
    // Make the A matrix Hermitian.
    // It sets Aji = conjf( Aij ) for j < i, that is, copy lower
    // triangle to upper triangle.
    //================================================================
    for (int i = 0; i < n; ++i) {
        A(i,i) = creal(A(i,i));
        for (int j = 0; j < i; ++j) {
            A(j,i) = conjf(A(i,j));
        }
    }
    int zerocol = param[PARAM_ZEROCOL].i;
    if (zerocol >= 0 && zerocol < n) {
        LAPACKE_claset_work(
            LAPACK_COL_MAJOR, 'F', n, 1, 0.0, 0.0, &A(0, zerocol), lda);
        LAPACKE_claset_work(
            LAPACK_COL_MAJOR, 'F', 1, n, 0.0, 0.0, &A(zerocol, 0), lda);
    }

    plasma_complex32_t *Aref = NULL;
    if (test) {
        Aref = (plasma_complex32_t*)malloc(
            (size_t)lda*n*sizeof(plasma_complex32_t));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(plasma_complex32_t));
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    int plainfo = plasma_chetrf(uplo, n, A, lda, ipiv, T, ldt, ipiv2);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_cpotrf(n) / time / 1e9;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        if (plainfo == 0) {
            // compute the residual norm ||A-bx||
            int nrhs = param[PARAM_NRHS].i;
            int ldb = imax(1, n + param[PARAM_PADB].i);
            int ldx = ldb;

            plasma_complex32_t *B = (plasma_complex32_t*)malloc(
                (size_t)ldb*nrhs*sizeof(plasma_complex32_t));
            assert(B != NULL);
            plasma_complex32_t *X = (plasma_complex32_t*)malloc(
                (size_t)ldx*nrhs*sizeof(plasma_complex32_t));
            assert(X != NULL);

            // set up right-hand-side B = A*rand
            retval = LAPACKE_clarnv(1, seed, (size_t)ldb*nrhs, X);
            assert(retval == 0);
            plasma_cgemm(PlasmaNoTrans, PlasmaNoTrans,
                         n, nrhs, n,
                         1.0, Aref, lda,
                                 X, ldx,
                         0.0,    B, ldb);

            // copy B to X
            LAPACKE_clacpy_work(
                LAPACK_COL_MAJOR, 'F', n, nrhs, B, ldb, X, ldx);

            // solve for X
            int iinfo = plasma_chetrs(
                uplo, n, nrhs, A, lda, ipiv, T, ldt, ipiv2, X, ldx);
            if (iinfo != 0) printf( " chetrs failed, info = %d\n", iinfo );

            // compute residual vector
            plasma_cgemm(PlasmaNoTrans, PlasmaNoTrans,
                         n, nrhs, n,
                         -1.0, Aref, lda,
                                  X, ldx,
                          1.0,    B, ldb);

            // compute various norms
            float *work = NULL;
            work = (float*)malloc((size_t)n*sizeof(float));
            assert(work != NULL);

            float Anorm = LAPACKE_clange_work(
                LAPACK_COL_MAJOR, 'F', n, n,    Aref, lda, work);
            float Xnorm = LAPACKE_clange_work(
                LAPACK_COL_MAJOR, 'I', n, nrhs, X, ldb, work);
            float Rnorm = LAPACKE_clange_work(
                LAPACK_COL_MAJOR, 'I', n, nrhs, B, ldb, work);
            float residual = Rnorm/(n*Anorm*Xnorm);

            param[PARAM_ERROR].d = residual;
            param[PARAM_SUCCESS].i = residual < tol;

            // free workspaces
            free(work);
            free(X);
            free(B);
        }
        else {
            // check for zero column, by updating zerocol with ipiv
            zerocol ++; // zerocol is 0-based, info is 1-based
            for (int i = nb; i < plainfo; i++) {
                if (ipiv[i] == zerocol) {
                    zerocol = i+1;
                }
                else if (i == zerocol) {
                    zerocol = ipiv[i];
                }
            }
            if (plainfo == zerocol) {
                param[PARAM_ERROR].d = 0.0;
                param[PARAM_SUCCESS].i = 1;
            }
            else {
                param[PARAM_ERROR].d = INFINITY;
                param[PARAM_SUCCESS].i = 0;
            }
        }
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A); free(ipiv); free(ipiv2);
    free(T);
    if (test)
        free(Aref);
}
