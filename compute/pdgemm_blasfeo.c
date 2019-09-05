/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @generated from compute/pzgemm.c, normal z -> d, Thu Aug  8 10:18:20 2019
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



//void plasma_core_omp_dgemm_blasfeo(
//    plasma_enum_t transa, plasma_enum_t transb,
//    int m, int n, int k,
//    double alpha, const double *A, int lda,
//                              const double *B, int ldb,
//    double beta,        double *C, int ldc,
//    double beta,        struct blasfeo_dmat *sC, int ci, int cj,
//    plasma_sequence_t *sequence, plasma_request_t *request);



#define A(m, n) (double*)plasma_tile_addr(A, m, n)
#define B(m, n) (double*)plasma_tile_addr(B, m, n)
#define C(m, n) (double*)plasma_tile_addr(C, m, n)

/***************************************************************************//**
 * Parallel tile matrix-matrix multiplication.
 * @see plasma_omp_dgemm
 ******************************************************************************/
void plasma_pdgemm_blasfeo(plasma_enum_t transa, plasma_enum_t transb,
                   double alpha, plasma_desc_t A,
                                             plasma_desc_t B,
                   double beta,  plasma_desc_t C,
                   plasma_sequence_t *sequence, plasma_request_t *request)
{
    // Return if failed sequence.
    if (sequence->status != PlasmaSuccess)
        return;
    if (A.type == PlasmaGeneral) {
        for (int m = 0; m < C.mt; m++) {
            int mvcm = plasma_tile_mview(C, m);
            for (int n = 0; n < C.nt; n++) {
#if HAVE_BLASFEO_API
				int sdcn = plasma_tile_nmain(C, n);
#else
				int ldcm = plasma_tile_mmain(C, m);
#endif
                int nvcn = plasma_tile_nview(C, n);
                //=========================================
                // alpha*A*B does not contribute; scale C
                //=========================================
                int inner_k = transa == PlasmaNoTrans ? A.n : A.m;
                if (alpha == 0.0 || inner_k == 0) {
                    int ldam = imax(1, plasma_tile_mmain(A, 0));
                    int ldbk = imax(1, plasma_tile_mmain(B, 0));
#if 0
                    plasma_core_omp_dgemm_blasfeo(
                        transa, transb,
                        mvcm, nvcn, 0,
                        alpha, A(0, 0), ldam,
                        B(0, 0), ldbk,
                        beta,  C(m, n), ldcm,
                        sequence, request);
#endif
                }
                else if (transa == PlasmaNoTrans) {
#if HAVE_BLASFEO_API
#else
					int ldam = plasma_tile_mmain(A, m);
#endif
                    //================================
                    // PlasmaNoTrans / PlasmaNoTrans
                    //================================
                    if (transb == PlasmaNoTrans) {
                        // printf("a:n b:n\n");
                        for (int k = 0; k < A.nt; k++) {
                            int nvak = plasma_tile_nview(A, k);
							int sdak = plasma_tile_nmain(A, k);
                            int sdbn = plasma_tile_nmain(B, n);
                            double zbeta = k == 0 ? beta : 1.0;
							// create blasfeo matrices
							struct blasfeo_dmat sA, sB, sC;
							blasfeo_create_dmat(m, k, &sA, A(m, k));
							sA.cn = sdak;
							blasfeo_create_dmat(k, n, &sB, B(k, n));
							sB.cn = sdbn;
							blasfeo_create_dmat(m, n, &sC, C(m, n));
							sC.cn = sdcn;
							// call LA routine
                            plasma_core_omp_dgemm_blasfeo(
                                transa, transb,
                                mvcm, nvcn, nvak,
                                alpha, &sA, 0, 0,
                                &sB, 0, 0,
                                zbeta, &sC, 0, 0,
                                sequence, request);
                        }
                    }
                    //=====================================
                    // PlasmaNoTrans / Plasma[_Conj]Trans
                    //=====================================
                    else {
                        // printf("a:n b:t\n");
                        for (int k = 0; k < A.nt; k++) {
                            int nvak = plasma_tile_nview(A, k);
							int sdak = plasma_tile_nmain(A, k);
                            int sdbn = plasma_tile_mmain(B, k); 
                            double zbeta = k == 0 ? beta : 1.0;
                            struct blasfeo_dmat sA, sB, sC;
							blasfeo_create_dmat(m, k, &sA, A(m, k));
							sA.cn = sdak;
							blasfeo_create_dmat(n, k, &sB, B(n, k));
							sB.cn = sdbn;
							blasfeo_create_dmat(m, n, &sC, C(m, n));
							sC.cn = sdcn;
                            plasma_core_omp_dgemm_blasfeo(
                                transa, transb,
                                mvcm, nvcn, nvak,
                                alpha, &sA, 0, 0,
                                &sB, 0, 0,
                                zbeta, &sC, 0, 0,
                                sequence, request);
                        }
                    }
                }
                //=====================================
                // Plasma[_Conj]Trans / PlasmaNoTrans
                //=====================================
                else {
                    // printf("a:t b:n\n");
                    if (transb == PlasmaNoTrans) {
                        for (int k = 0; k < A.mt; k++) {

							int sdak = plasma_tile_mmain(A, m);
                            int sdbn = plasma_tile_nmain(B, n);
                            int mvak = plasma_tile_mview(A, k);
                            double zbeta = k == 0 ? beta : 1.0;
                            struct blasfeo_dmat sA, sB, sC;
							blasfeo_create_dmat(k, m, &sA, A(k, m));
							sA.cn = sdak;
							blasfeo_create_dmat(k, n, &sB, B(k, n));
							sB.cn = sdbn;
							blasfeo_create_dmat(m, n, &sC, C(m, n));
							sC.cn = sdcn;
                            plasma_core_omp_dgemm_blasfeo(
                                transa, transb,
                                mvcm, nvcn, mvak,
                                alpha, &sA, 0, 0,
                                &sB, 0, 0,
                                zbeta, &sC, 0, 0,
                                sequence, request);
                        }
                    }
                    //==========================================
                    // Plasma[_Conj]Trans / Plasma[_Conj]Trans
                    //==========================================
                    else {
                        // printf("a:t b:t\n");
                        for (int k = 0; k < A.mt; k++) {

                            int mvak = plasma_tile_mview(A, k);
							int sdak = plasma_tile_mmain(A, m);
                            int sdbn = plasma_tile_mmain(B, k);
                            double zbeta = k == 0 ? beta : 1.0;
                            struct blasfeo_dmat sA, sB, sC;
							blasfeo_create_dmat(k, m, &sA, A(k, m));
							sA.cn = sdak;
							blasfeo_create_dmat(n, k, &sB, B(n, k));
							sB.cn = sdbn;
							blasfeo_create_dmat(m, n, &sC, C(m, n));
							sC.cn = sdcn;
                            plasma_core_omp_dgemm_blasfeo(
                                transa, transb,
                                mvcm, nvcn, mvak,
                                alpha, &sA, 0, 0,
                                &sB, 0, 0,
                                zbeta, &sC, 0, 0,
                                sequence, request);
                        }
                    }
                }
            }
        }
    }
    else if (A.type == PlasmaGeneralBand) {
        for (int m = 0; m < C.mt; m++) {
            int mvcm = plasma_tile_mview(C, m);
            int ldcm = plasma_tile_mmain(C, m);
            for (int n = 0; n < C.nt; n++) {
                int nvcn = plasma_tile_nview(C, n);
                //=========================================
                // alpha*A*B does not contribute; scale C
                //=========================================
                int inner_k = transa == PlasmaNoTrans ? A.n : A.m;
                if (alpha == 0.0 || inner_k == 0) {
                    int ldam = imax(1, plasma_tile_mmain(A, 0));
                    int ldbk = imax(1, plasma_tile_mmain(B, 0));
#if 0
                    plasma_core_omp_dgemm_blasfeo(
                        transa, transb,
                        mvcm, nvcn, 0,
                        alpha, A(0, 0), ldam,
                        B(0, 0), ldbk,
                        beta,  C(m, n), ldcm,
                        sequence, request);
#endif
                }
                else if (transa == PlasmaNoTrans) {
                    //================================
                    // PlasmaNoTrans / PlasmaNoTrans
                    //================================
                    if (transb == PlasmaNoTrans) {
                        int k_start = (imax(0, m*A.mb-A.kl)) / A.nb;
                        int k_end = (imin(A.n-1, (m+1)*A.mb+A.ku-1)) / A.nb;
                        //printf("[%s]: m=%d\tn=%d\tk_s=%d\tk_e=%d\n",
                        //       __FILE__, m, n, k_start, k_end);
                        for (int k = k_start; k <= k_end; k++) {
                            int ldam = plasma_tile_mmain_band(A, m, k);
                            int nvak = plasma_tile_nview(A, k);
                            int ldbk = plasma_tile_mmain(B, k);
                            double zbeta = k == 0 ? beta : 1.0;

#if 0
                            plasma_core_omp_dgemm_blasfeo(
                                transa, transb,
                                mvcm, nvcn, nvak,
                                alpha, A(m, k), ldam,
                                B(k, n), ldbk,
                                zbeta, C(m, n), ldcm,
                                sequence, request);
#endif
                        }
                    }
                    //=====================================
                    // PlasmaNoTrans / Plasma[_Conj]Trans
                    //=====================================
                    else {
                        assert(0); // not implemented
                    }
                }
                else {
                    //=====================================
                    // Plasma[_Conj]Trans / PlasmaNoTrans
                    //=====================================
                    if (transb == PlasmaNoTrans) {
                        assert(0); // not implemented
                    }
                    //==========================================
                    // Plasma[_Conj]Trans / Plasma[_Conj]Trans
                    //==========================================
                    else {
                        assert(0); // not implemented
                    }
                }
            }
        }
    }
}
