/**
 *
 * @file pzgeqrf.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley, Univ. of Colorado Denver and
 *  Univ. of Manchester.
 *
 * @version 3.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Sistek
 * @date 2016-7-6
 * @precisions normal z -> s d c
 *
 **/

#include "plasma_async.h"
#include "plasma_descriptor.h"
#include "plasma_types.h"
#include "plasma_internal.h"
#include "core_blas_z.h"

#define A(m, n) ((PLASMA_Complex64_t*) plasma_getaddr(A, m, n))
#define T(m, n) ((PLASMA_Complex64_t*) plasma_getaddr(T, m, n))
/***************************************************************************//**
 *  Parallel tile QR factorization - dynamic scheduling
 * @see PLASMA_zgeqrf_Tile_Async
 **/
void plasma_pzgeqrf(PLASMA_desc A, PLASMA_desc T,
                    PLASMA_sequence *sequence, PLASMA_request *request)
{
    int k, m, n;
    int ldak, ldam;
    int tempkm, tempkn, tempnn, tempmm;
    int ib;

    if (sequence->status != PLASMA_SUCCESS)
        return;

    // Jakub S.: How to set inner blocking?
    ib = A.nb;

    for (k = 0; k < imin(A.mt, A.nt); k++) {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        ldak = BLKLDD(A, k);
        CORE_OMP_zgeqrt(
            tempkm, tempkn, ib, T.nb,
            A(k, k), ldak,
            T(k, k), T.mb);

        for (n = k+1; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            // Plasma_ConjTrans will be converted to PlasmaTrans in
            // automatic datatype conversion, which is what we
            // want here.
            // PlasmaConjTrans is protected from this conversion.
            CORE_OMP_zunmqr(
                PlasmaLeft, Plasma_ConjTrans,
                tempkm, tempnn, tempkm, ib, T.nb,
                A(k, k), ldak,
                T(k, k), T.mb,
                A(k, n), ldak);
        }
        for (m = k+1; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            CORE_OMP_ztsqrt(
                tempmm, tempkn, ib, T.nb,
                A(k, k), ldak,
                A(m, k), ldam,
                T(m, k), T.mb);

            for (n = k+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                CORE_OMP_ztsmqr(
                    PlasmaLeft, Plasma_ConjTrans,
                    A.mb, tempnn, tempmm, tempnn, A.nb, ib, T.nb,
                    A(k, n), ldak,
                    A(m, n), ldam,
                    A(m, k), ldam,
                    T(m, k), T.mb);
            }
        }
    }
}
