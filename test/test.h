/**
 *
 * @file test.h
 *
 *  PLASMA testing harness.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver.
 *
 * @version 3.0.0
 * @author Jakub Kurzak
 * @date 2016-01-01
 *
 **/
#ifndef TEST_H
#define TEST_H

//==============================================================================
// parameter labels
//==============================================================================
typedef enum {
    //------------------------------------------------------
    // input parameters
    //------------------------------------------------------
    PARAM_ITER,   // outer product iteration?
    PARAM_OUTER,  // outer product iteration?
    PARAM_TEST,   // test the solution?
    PARAM_TOL,    // tolerance
    PARAM_TRANSA, // transposition of A
    PARAM_TRANSB, // transposition of B
    PARAM_M,      // M dimension
    PARAM_N,      // N dimension
    PARAM_K,      // K dimension
    PARAM_PADA,   // padding of A
    PARAM_PADB,   // padding of B
    PARAM_PADC,   // padding of C

    //------------------------------------------------------
    // output parameters
    //------------------------------------------------------
    PARAM_SUCCESS, // success indicator
    PARAM_ERROR,   // numerical error
    PARAM_TIME,    // time to solution
    PARAM_GFLOPS,  // GFLOPS rate

    //------------------------------------------------------
    // Keep at the end!
    //------------------------------------------------------
    PARAM_SIZEOF   // size of parameter array
} param_label_t;

//==============================================================================
// parameter descriptions
//==============================================================================
static const char *ParamUsage[][2] = {
    {"--iter=", "number of iterations per set of parameters [default: 1]"},
    {"--outer=[y|n]", "outer product iteration [default: n]"},
    {"--test=[y|n]", "test the solution [default: y]"},
    {"--tol=", "tolerance [default: 50]"},
    {"--transa=[n|t|c]", "transposition of A [default: n]"},
    {"--transb=[n|t|c]", "transposition of B [default: n]"},
    {"--m=", "M dimension (number of rows) [default: 1000]"},
    {"--n=", "N dimension (number of columns) [default: 1000]"},
    {"--k=", "K dimension (number of rows or columns) [default: 1000]"},
    {"--pada=", "padding of A [default: 0]"},
    {"--padb=", "padding of B [default: 0]"},
    {"--padc=", "padding of C [default: 0]"}
};

//==============================================================================
// tester infrastructure
//==============================================================================
// parameter value type
typedef union {
    int i;          // integer
    char c;         // character
    double d;       // double precision
} param_value_t;

// parameter type
typedef struct {
    int num;            // number of values for a parameter
    int pos;            // current position in the array
    int size;           // size of parameter values array
    param_value_t *val; // array of values for a parameter
} param_t;

// hiding double from precision translation when used for taking time
typedef double plasma_time_t;

// initial size of values array
static const int InitValArraySize = 1024;

// indentation of option descriptions
static const int DescriptionIndent = -20;

// maximum length of info string
static const int InfoLen = 1024;

// spacing in info output string
static const int InfoSpacing = 12;

// function declarations
void print_main_usage();
void print_routine_usage(char *name);
void print_usage(int label);
void test_routine(char *name, param_value_t param[]);
void time_routine(char *name, param_value_t pval[]);
void run_routine(char *name, param_value_t pval[], char *info);
void param_init(param_t param[]);
int param_read(int argc, char **argv, param_t param[]);
int param_starts_with(char *str, char *prefix);
void param_scan_int(char *str, param_t *param);
void param_scan_char(char *str, param_t *param);
void param_scan_double(char *str, param_t *param);
void param_add_int(int val, param_t *param);
void param_add_char(char cval, param_t *param);
void param_add_double(double dval, param_t *param);
int param_step_inner(param_t param[]);
int param_step_outer(param_t param[], int idx);
int param_snap(param_t param[], param_value_t value[]);

#include "test_s.h"
#include "test_d.h"
#include "test_c.h"
#include "test_z.h"

#endif // TEST_H