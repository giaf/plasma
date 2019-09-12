/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzpotrf.c, normal z -> d, Thu Aug  8 17:08:17 2019
 *
 **/

#include "plasma_async.h"
#include "plasma_context.h"
#include "plasma_descriptor.h"
#include "plasma_internal.h"
#include "plasma_types.h"
#include "plasma_workspace.h"
#include <plasma_core_blas.h>

#include "blasfeo_d_aux.h"

#define A(m, n) (double*)plasma_tile_addr(A, m, n)

/***************************************************************************//**
 *  Parallel tile Cholesky factorization.
 * @see plasma_omp_dpotrf
 ******************************************************************************/
void plasma_pdpotrf_blasfeo(plasma_enum_t uplo, plasma_desc_t A,
                    plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;

    //==============
    // PlasmaLower
    //==============
    if (uplo == PlasmaLower) {
        fprintf(stderr,"before malloc\n");
//        double *ptr = malloc(A.gmt*A.mb*sizeof(double)); //do we need to use dA for other opertations too?
		double *ptr = A.vector;
        
        for (int k = 0; k < A.mt; k++)
		{
            int mvak = plasma_tile_mview(A, k);
            int ldak = plasma_tile_mmain(A, k);
            int sdak = plasma_tile_nmain(A, k);

            // create blasfeo matrix
            struct blasfeo_dmat sA1;
            blasfeo_create_dmat(ldak, sdak, &sA1, A(k,k));
            sA1.dA = ptr + k*A.nb; // TODO map vector 1-to-1 to the single tiles in the same order
			sA1.use_dA = mvak;
#if 1
            plasma_core_omp_dpotrf_blasfeo(
                PlasmaLower, mvak,
                &sA1, 0, 0,
                A.nb*k,
                sequence, request);
#endif
//goto end;

            for (int m = k+1; m < A.mt; m++)
			{
                int mvam = plasma_tile_mview(A, m);
                int ldam = plasma_tile_mmain(A, m);

                struct blasfeo_dmat sA2;
                fprintf(stderr, "before create dmat1\n");
                blasfeo_create_dmat(ldam, sdak, &sA2, A(m,k));

#if 1
                plasma_core_omp_dtrsm_blasfeo(
                    PlasmaRight, PlasmaLower,
                    PlasmaConjTrans, PlasmaNonUnit,
                    mvam, A.mb,
                    1.0, &sA1, 0, 0,
                         &sA2, 0, 0,
                    sequence, request);
#endif
            }
            for (int m = k+1; m < A.mt; m++) {
                int mvam = plasma_tile_mview(A, m);
                int ldam = plasma_tile_mmain(A, m);
                int sdam = plasma_tile_nmain(A, m);

                struct blasfeo_dmat sA2, sA3;
                blasfeo_create_dmat(ldam, sdak, &sA2, A(m,k));
                blasfeo_create_dmat(ldam, sdam, &sA3, A(m,m));

                printf("before dsyrk\n");
#if 1
                plasma_core_omp_dsyrk_blasfeo(
                    PlasmaLower, PlasmaNoTrans,
                    mvam, A.mb,
                    -1.0, &sA2, 0, 0,
                     1.0, &sA3, 0, 0,
                    sequence, request);
#endif

                for (int n = k+1; n < m; n++) {
                    int ldan = plasma_tile_mmain(A, n);
                    int sdan = plasma_tile_nmain(A, n);

                    struct blasfeo_dmat sA4, sA5;
                    blasfeo_create_dmat(ldan, sdak, &sA4, A(n,k));
                    blasfeo_create_dmat(mvam, sdan, &sA5, A(m,n));

#if 1
                    plasma_core_omp_dgemm_blasfeo(
                        PlasmaNoTrans, PlasmaConjTrans,
                        mvam, A.mb, A.mb,
                        -1.0, &sA2, 0, 0,
                              &sA4, 0, 0,
                         1.0, &sA5, 0, 0,
                        sequence, request);
#endif
                }
            }
        }
//end:
//        free(ptr);
    }
    //==============
    // PlasmaUpper
    //==============
    #if 0
    else {
        for (int k = 0; k < A.nt; k++) {
            int nvak = plasma_tile_nview(A, k);
            int ldak = plasma_tile_mmain(A, k);
            plasma_core_omp_dpotrf_blasfeo(
                PlasmaUpper, nvak,
                A(k, k), ldak,
                A.nb*k,
                sequence, request);

            for (int m = k+1; m < A.nt; m++) {
                int nvam = plasma_tile_nview(A, m);
                plasma_core_omp_dtrsm(
                    PlasmaLeft, PlasmaUpper,
                    PlasmaConjTrans, PlasmaNonUnit,
                    A.nb, nvam,
                    1.0, A(k, k), ldak,
                         A(k, m), ldak,
                    sequence, request);
            }
            for (int m = k+1; m < A.nt; m++) {
                int nvam = plasma_tile_nview(A, m);
                int ldam = plasma_tile_mmain(A, m);
                plasma_core_omp_dsyrk(
                    PlasmaUpper, PlasmaConjTrans,
                    nvam, A.mb,
                    -1.0, A(k, m), ldak,
                     1.0, A(m, m), ldam,
                    sequence, request);

                for (int n = k+1; n < m; n++) {
                    int ldan = plasma_tile_mmain(A, n);
                    plasma_core_omp_dgemm(
                        PlasmaConjTrans, PlasmaNoTrans,
                        A.mb, nvam, A.mb,
                        -1.0, A(k, n), ldak,
                              A(k, m), ldak,
                         1.0, A(n, m), ldan,
                        sequence, request);
                }
            }
        }
    }
    #endif
}
