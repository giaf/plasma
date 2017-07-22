/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from test/test_zlangb.c, normal z -> s, Sat Jul 22 09:25:37 2017
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
float slangb_(char *, int *, int *, int *, float *, int *, float *);
 //   float value = plasma_slangb(norm, m, n, kl, ku,  AB, ldab);
float plasma_slangb(plasma_enum_t, int, int, int, int, float *, int);
/***************************************************************************//**
 *
 * @brief Tests SLANGB.
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets flags in param indicating which parameters are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_slangb(param_value_t param[], bool run)
{
    // pwu: experiements:
    //test_lapack_slangb();
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_NORM   ].used = true;
    param[PARAM_DIM    ].used = PARAM_USE_M | PARAM_USE_N;
    param[PARAM_PADA   ].used = true;
    param[PARAM_NB     ].used = true;
    param[PARAM_KL     ].used = true;
    param[PARAM_KU     ].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t norm = plasma_norm_const(param[PARAM_NORM].c);

    int m = param[PARAM_DIM].dim.m;
    int n = param[PARAM_DIM].dim.n;
    int kl = param[PARAM_KL].i;
    int ku = param[PARAM_KU].i;

    int lda = imax(1, m + param[PARAM_PADA].i);

    int test = param[PARAM_TEST].c == 'y';
    float eps = LAPACKE_slamch('E');

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    float *A =
        (float*)malloc((size_t)lda*n*sizeof(float));
    assert(A != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_slarnv(1, seed, (size_t)lda*n, A);
    assert(retval == 0);
    // zero out elements outside the band
    for (int i = 0; i < m; i++) {
        for (int j = i+ku+1; j < n; j++) A[i + j*lda] = 0.0;
    }
    for (int j = 0; j < n; j++) {
        for (int i = j+kl+1; i < m; i++) A[i + j*lda] = 0.0;
    }
    float *Aref = NULL;
    if (test) {
        Aref = (float*)malloc(
            (size_t)lda*n*sizeof(float));
        assert(Aref != NULL);

        memcpy(Aref, A, (size_t)lda*n*sizeof(float));
    }

    int nb = param[PARAM_NB].i;
    // band matrix A in skewed LAPACK storage
    int kut  = (ku+kl+nb-1)/nb; // # of tiles in upper band (not including diagonal)
    int klt  = (kl+nb-1)/nb;    // # of tiles in lower band (not including diagonal)
    int ldab = (kut+klt+1)*nb;  // since we use sgetrf on panel, we pivot back within panel.
                                // this could fill the last tile of the panel,
                                // and we need extra NB space on the bottom
    float *AB = NULL;
    AB = (float*)malloc((size_t)ldab*n*sizeof(float));
    assert(AB != NULL);
    // convert into LAPACK's skewed storage
    for (int j = 0; j < n; j++) {
        int i_kl = imax(0,   j-ku);
        int i_ku = imin(m-1, j+kl);
        for (int i = 0; i < ldab; i++) AB[i + j*ldab] = 0.0;
        for (int i = i_kl; i <= i_ku; i++) AB[kl + i-(j-ku) + j*ldab] = A[i + j*lda];
    }

    float *ABref = NULL;
    if (test) {
        ABref = (float*)malloc(
            (size_t)ldab*n*sizeof(float));
        assert(ABref != NULL);

        memcpy(ABref, AB, (size_t)ldab*n*sizeof(float));
    }
    
    //================================================================
    // Run and time PLASMA.
    //================================================================
    plasma_time_t start = omp_get_wtime();
    float value = plasma_slangb(norm, m, n, kl, ku,  AB, ldab);
    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;

    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_slange(m, n, norm) / time / 1e9;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        char *cnorm;
	float *workspace=NULL;
        switch (norm) {
        case PlasmaMaxNorm:
            cnorm = "m";
            break;
        case PlasmaOneNorm:
            cnorm = "1";
            break;
	case PlasmaInfNorm:
	    cnorm = "i";
	    workspace = (float*) malloc(n*sizeof(float));
	    break;
	case PlasmaFrobeniusNorm:
	    cnorm = "f";
	    break;
        default:
            assert(0);
        }

        float valueRef =
            //LAPACKE_slange(LAPACK_COL_MAJOR, lapack_const(norm),
            //               kl+ku+1, n, ABref+kl, ldab);
            //slange_(cnorm, &kuu, &n, ABref+kl, &ldab, NULL);
            slangb_(cnorm, &n,&kl, &ku, ABref+kl, &ldab, workspace);

        // Calculate relative error
        float error = fabsf(value-valueRef);

        if (valueRef != 0)
            error /= valueRef;
        float tol = eps;
        float normalize = 1;
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
    free(A);free(AB);
    if (test) {
        free(Aref);free(ABref);
    }
}



