/**
 *
 * @file test_z.h
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef TEST_Z_H
#define TEST_Z_H

//==============================================================================
// test routines
//==============================================================================
void test_zgemm(param_value_t param[], char *info);
void test_zhemm(param_value_t param[], char *info);
void test_zher2k(param_value_t param[], char *info);
void test_zherk(param_value_t param[], char *info);
void test_zposv(param_value_t param[], char *info);
void test_zpotrf(param_value_t param[], char *info);
void test_zpotrs(param_value_t param[], char *info);
void test_zsymm(param_value_t param[], char *info);
void test_zsyr2k(param_value_t param[], char *info);
void test_zsyrk(param_value_t param[], char *info);
void test_ztrsm(param_value_t param[], char *info);

#endif // TEST_Z_H
