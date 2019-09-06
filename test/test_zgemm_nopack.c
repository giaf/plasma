/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/
#include "test.h"
#include "flops.h"
#include "plasma.h"
#include "core_lapack.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <omp.h>

#define COMPLEX

/***************************************************************************//**
 *
 * @brief Tests ZGEMM.
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets used flags in param indicating parameters that are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_zgemm_nopack(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_TRANSA ].used = true;
    param[PARAM_TRANSB ].used = true;
    param[PARAM_DIM    ].used = PARAM_USE_M | PARAM_USE_N | PARAM_USE_K;
    param[PARAM_ALPHA  ].used = true;
    param[PARAM_BETA   ].used = true;
    param[PARAM_PADA   ].used = true;
    param[PARAM_PADB   ].used = true;
    param[PARAM_PADC   ].used = true;
    param[PARAM_NB     ].used = true;
    if (! run)
        return;

    //================================================================
    // Set parameters.
    //================================================================
    plasma_enum_t transa = plasma_trans_const(param[PARAM_TRANSA].c);
    plasma_enum_t transb = plasma_trans_const(param[PARAM_TRANSB].c);

    int m = param[PARAM_DIM].dim.m;
    int n = param[PARAM_DIM].dim.n;
    int k = param[PARAM_DIM].dim.k;

    int Am, An;
    int Bm, Bn;
    int Cm, Cn;

    if (transa == PlasmaNoTrans) {
        Am = m;
        An = k;
    }
    else {
        Am = k;
        An = m;
    }
    if (transb == PlasmaNoTrans) {
        Bm = k;
        Bn = n;
    }
    else {
        Bm = n;
        Bn = k;
    }
    Cm = m;
    Cn = n;

    int lda = imax(1, Am + param[PARAM_PADA].i);
    int ldb = imax(1, Bm + param[PARAM_PADB].i);
    int ldc = imax(1, Cm + param[PARAM_PADC].i);

    int test = param[PARAM_TEST].c == 'y';
    double eps = LAPACKE_dlamch('E');

#ifdef COMPLEX
    plasma_complex64_t alpha = param[PARAM_ALPHA].z;
    plasma_complex64_t beta  = param[PARAM_BETA].z;
#else
    double alpha = creal(param[PARAM_ALPHA].z);
    double beta  = creal(param[PARAM_BETA].z);
#endif

    //================================================================
    // Set tuning parameters.
    //================================================================
    plasma_set(PlasmaTuning, PlasmaDisabled);
    plasma_set(PlasmaNb, param[PARAM_NB].i);

    //================================================================
    // Allocate and initialize arrays.
    //================================================================
    plasma_complex64_t *A =
        (plasma_complex64_t*)malloc((size_t)lda*An*sizeof(plasma_complex64_t));
    assert(A != NULL);

    plasma_complex64_t *B =
        (plasma_complex64_t*)malloc((size_t)ldb*Bn*sizeof(plasma_complex64_t));
    assert(B != NULL);

    plasma_complex64_t *C =
        (plasma_complex64_t*)malloc((size_t)ldc*Cn*sizeof(plasma_complex64_t));
    assert(C != NULL);

    int seed[] = {0, 0, 0, 1};
    lapack_int retval;
    retval = LAPACKE_zlarnv(1, seed, (size_t)lda*An, A);
    assert(retval == 0);

    retval = LAPACKE_zlarnv(1, seed, (size_t)ldb*Bn, B);
    assert(retval == 0);

    retval = LAPACKE_zlarnv(1, seed, (size_t)ldc*Cn, C);
    assert(retval == 0);

    plasma_complex64_t *Cref = NULL;
    if (test) {
        Cref = (plasma_complex64_t*)malloc(
            (size_t)ldc*Cn*sizeof(plasma_complex64_t));
        assert(Cref != NULL);

        memcpy(Cref, C, (size_t)ldc*Cn*sizeof(plasma_complex64_t));
    }

    //================================================================
    // Run and time PLASMA.
    //================================================================

/*     plasma_zgemm(
        transa, transb,
        m, n, k,
        alpha, A, lda,
               B, ldb,
         beta, C, ldc);
 */   
    // Get PLASMA context.
    plasma_context_t *plasma = plasma_context_self();
    if (plasma == NULL) {
        plasma_error("PLASMA not initialized");
        return PlasmaErrorNotInitialized;
    }

    // Check input arguments.
    if ((transa != PlasmaNoTrans) &&
        (transa != PlasmaTrans) &&
        (transa != PlasmaConjTrans)) {
        plasma_error("illegal value of transa");
        return -1;
    }
    if ((transb != PlasmaNoTrans) &&
        (transb != PlasmaTrans) &&
        (transb != PlasmaConjTrans)) {
        plasma_error("illegal value of transb");
        return -2;
    }
    if (m < 0) {
        plasma_error("illegal value of m");
        return -3;
    }
    if (n < 0) {
        plasma_error("illegal value of n");
        return -4;
    }
    if (k < 0) {
        plasma_error("illegal value of k");
        return -5;
    }

    int am, an;
    int bm, bn;
    if (transa == PlasmaNoTrans) {
        am = m;
        an = k;
    }
    else {
        am = k;
        an = m;
    }
    if (transb == PlasmaNoTrans) {
        bm = k;
        bn = n;
    }
    else {
        bm = n;
        bn = k;
    }

    if (lda < imax(1, am)) {
        plasma_error("illegal value of lda");
        return -8;
    }
    if (ldb < imax(1, bm)) {
        plasma_error("illegal value of ldb");
        return -10;
    }
    if (ldc < imax(1, m)) {
        plasma_error("illegal value of ldc");
        return -13;
    }

    // quick return
    if (m == 0 || n == 0 || ((alpha == 0.0 || k == 0) && beta == 1.0))
        return PlasmaSuccess;

    // Tune parameters.
    if (plasma->tuning)
        plasma_tune_gemm(plasma, PlasmaComplexDouble, m, n, k);

    // Set tiling parameters.
    int nb = plasma->nb;

    // Create tile matrices.
    plasma_desc_t AA;
    plasma_desc_t BB;
    plasma_desc_t CC;
    int retval1;
    retval1 = plasma_desc_general_create(PlasmaComplexDouble, nb, nb,
                                        am, an, 0, 0, am, an, &AA);
    if (retval1 != PlasmaSuccess) {
        plasma_error("plasma_desc_general_create() failed");
        return retval1;
    }
    retval1 = plasma_desc_general_create(PlasmaComplexDouble, nb, nb,
                                        bm, bn, 0, 0, bm, bn, &BB);
    if (retval1 != PlasmaSuccess) {
        plasma_error("plasma_desc_general_create() failed");
        plasma_desc_destroy(&AA);
        return retval1;
    }
    retval1 = plasma_desc_general_create(PlasmaComplexDouble, nb, nb,
                                        m, n, 0, 0, m, n, &CC);
    if (retval1 != PlasmaSuccess) {
        plasma_error("plasma_desc_general_create() failed");
        plasma_desc_destroy(&AA);
        plasma_desc_destroy(&BB);
        return retval1;
    }

    // Initialize sequence.
    plasma_sequence_t sequence;
    retval1 = plasma_sequence_init(&sequence);

    // Initialize request.
    plasma_request_t request;
    retval1 = plasma_request_init(&request);

    // asynchronous block
    #pragma omp parallel
    #pragma omp master
    {
        // Translate to tile layout.
        plasma_omp_zge2desc(A, lda, AA, &sequence, &request);
        plasma_omp_zge2desc(B, ldb, BB, &sequence, &request);
        plasma_omp_zge2desc(C, ldc, CC, &sequence, &request);
    }
    
    plasma_time_t start = omp_get_wtime();
    #pragma omp parallel
    #pragma omp master
    {
        // Call the tile async function.
        plasma_omp_zgemm(transa, transb,
                         alpha, AA,
                                BB,
                         beta,  CC,
                         &sequence, &request);
        // Translate back to LAPACK layout.
        plasma_omp_zdesc2ge(CC, C, ldc, &sequence, &request);
    }
    // implicit synchronization

    plasma_time_t stop = omp_get_wtime();
    plasma_time_t time = stop-start;
    param[PARAM_TIME].d = time;
    param[PARAM_GFLOPS].d = flops_zgemm(m, n, k) / time / 1e9;
    
    #pragma omp parallel
    #pragma omp master
    {
    // Translate back to LAPACK layout.
        plasma_omp_zdesc2ge(CC, C, ldc, &sequence, &request);
    }
    // Free matrices in tile layout.
    plasma_desc_destroy(&AA);
    plasma_desc_destroy(&BB);
    plasma_desc_destroy(&CC);

    
    //plasma_time_t time = stop-start;

    //param[PARAM_TIME].d = time;
    //param[PARAM_GFLOPS].d = flops_zgemm(m, n, k) / time / 1e9;

    //================================================================
    // Test results by comparing to a reference implementation.
    //================================================================
    if (test) {
        // |R - R_ref|_p < gamma_{k+2} * |alpha| * |A|_p * |B|_p +
        //                 gamma_2 * |beta| * |C|_p
        // holds component-wise or with |.|_p as 1, inf, or Frobenius norm.
        // gamma_k = k*eps / (1 - k*eps), but we use
        // gamma_k = sqrt(k)*eps as a statistical average case.
        // Using 3*eps covers complex arithmetic.
        // See Higham, Accuracy and Stability of Numerical Algorithms, ch 2-3.
        double work[1];
        double Anorm = LAPACKE_zlange_work(
                           LAPACK_COL_MAJOR, 'F', Am, An, A,    lda, work);
        double Bnorm = LAPACKE_zlange_work(
                           LAPACK_COL_MAJOR, 'F', Bm, Bn, B,    ldb, work);
        double Cnorm = LAPACKE_zlange_work(
                           LAPACK_COL_MAJOR, 'F', Cm, Cn, Cref, ldc, work);

        cblas_zgemm(
            CblasColMajor,
            (CBLAS_TRANSPOSE)transa, (CBLAS_TRANSPOSE)transb,
            m, n, k,
            CBLAS_SADDR(alpha), A, lda,
                                B, ldb,
             CBLAS_SADDR(beta), Cref, ldc);

        plasma_complex64_t zmone = -1.0;
        cblas_zaxpy((size_t)ldc*Cn, CBLAS_SADDR(zmone), Cref, 1, C, 1);

        double error = LAPACKE_zlange_work(
                           LAPACK_COL_MAJOR, 'F', Cm, Cn, C,    ldc, work);
        double normalize = sqrt((double)k+2) * cabs(alpha) * Anorm * Bnorm
                         + 2 * cabs(beta) * Cnorm;
        if (normalize != 0)
            error /= normalize;

        param[PARAM_ERROR].d = error;
        param[PARAM_SUCCESS].i = error < 3*eps;
    }

    //================================================================
    // Free arrays.
    //================================================================
    free(A);
    free(B);
    free(C);
    if (test)
        free(Cref);
}
