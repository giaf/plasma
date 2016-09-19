/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef ICL_PLASMA_TYPES_H
#define ICL_PLASMA_TYPES_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
#ifdef PLASMA_WITH_MKL
#define lapack_complex_float plasma_complex32_t
#define lapack_complex_double plasma_complex64_t
#endif

/******************************************************************************/
#ifndef CBLAS_SADDR
#define CBLAS_SADDR(var) &(var)
#endif

/******************************************************************************/
enum {
    PlasmaByte          = 0,
    PlasmaInteger       = 1,
    PlasmaRealFloat     = 2,
    PlasmaRealDouble    = 3,
    PlasmaComplexFloat  = 4,
    PlasmaComplexDouble = 5
};

/***************************************************************************//**
 *
 *  PLASMA constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 **/
enum {
    PlasmaInvalid    = -1,

    PlasmaNoTrans    = 111,
    PlasmaTrans      = 112,
    PlasmaConjTrans  = 113,
    Plasma_ConjTrans = PlasmaConjTrans,

    PlasmaUpper      = 121,
    PlasmaLower      = 122,
    PlasmaFull       = 123,  // formerly PlasmaUpperLower

    PlasmaNonUnit    = 131,
    PlasmaUnit       = 132,

    PlasmaLeft       = 141,
    PlasmaRight      = 142,

    PlasmaForward    = 391,
    PlasmaBackward   = 392,

    PlasmaColumnwise = 401,
    PlasmaRowwise    = 402,

    PlasmaW          = 501,
    PlasmaA2         = 502
};

enum {
    PlasmaSuccess = 0,
    PlasmaErrorNotInitialized,
    PlasmaErrorNotSupported,
    PlasmaErrorIllegalValue,
    PlasmaErrorOutOfMemory,
    PlasmaErrorNullParameter,
    PlasmaErrorInternal,
    PlasmaErrorSequence
};

enum {
    PlasmaInplace,
    PlasmaOutplace
};

enum {
    PlasmaNb,
    PlasmaIb,
    PlasmaInplaceOutplace
};

/******************************************************************************/
typedef int plasma_enum_t;

typedef float  _Complex plasma_complex32_t;
typedef double _Complex plasma_complex64_t;

/******************************************************************************/
plasma_enum_t plasma_trans_const_t(char lapack_char);
plasma_enum_t plasma_uplo_const_t(char lapack_char);
plasma_enum_t plasma_diag_const_t(char lapack_char);
plasma_enum_t plasma_side_cons_t(char lapack_char);
plasma_enum_t plasma_direct_const_t(char lapack_char);
plasma_enum_t plasma_storev_const_t(char lapack_char);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_PLASMA_TYPES_H
