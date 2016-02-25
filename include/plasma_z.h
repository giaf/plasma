/**
 *
 * @file plasma_z.h
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver.
 *
 * @version 3.0.0
 * @author Jakub Kurzak
 * @date 2016-01-01
 * @precisions normal z -> s d c
 *
 **/
#ifndef PLASMA_Z_H
#define PLASMA_Z_H

/***************************************************************************//**
 *  Standard interface.
 **/
int PLASMA_zgemm(PLASMA_enum transA, PLASMA_enum transB,
                 int m, int n, int k,
                 PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int lda,
                                           PLASMA_Complex64_t *B, int ldb,
                 PLASMA_Complex64_t beta,  PLASMA_Complex64_t *C, int ldc);

/***************************************************************************//**
 *  Tile interface.
 **/
int PLASMA_zgemm_Tile(PLASMA_enum transA, PLASMA_enum transB,
                      PLASMA_Complex64_t alpha, PLASMA_desc *A,
                                                PLASMA_desc *B,
                      PLASMA_Complex64_t beta,  PLASMA_desc *C);

/***************************************************************************//**
 *  Tile asynchronous interface.
 **/
int PLASMA_zgemm_Tile_Async(PLASMA_enum transA, PLASMA_enum transB,
                            PLASMA_Complex64_t alpha, PLASMA_desc *A,
                                                      PLASMA_desc *B,
                            PLASMA_Complex64_t beta,  PLASMA_desc *C,
                            PLASMA_sequence *sequence, PLASMA_request *request);

#endif // PLASMA_Z_H
