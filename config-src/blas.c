
#if defined( MKL_ILP64 ) || defined( ILP64 )
    typedef long long myint;
#else
    typedef int myint;
#endif

#if defined( LOWERCASE )
    #define FORTRAN_NAME( lower, UPPER ) lower
#elif defined( UPPERCASE )       
    #define FORTRAN_NAME( lower, UPPER ) UPPER
#else                            
    #define FORTRAN_NAME( lower, UPPER ) lower ## _
#endif

#ifdef __cplusplus
    #define EXTERN_C extern "C"
#else
    #define EXTERN_C
#endif

#include <stdio.h>

#define sgemm FORTRAN_NAME( sgemm, SGEMM )

EXTERN_C
void sgemm( const char* transA, const char* transB,
            const myint* m, const myint* n, const myint* k,
            const float* alpha,
            const float* A, const myint* lda,
            const float* B, const myint* ldb,
            const float* beta,
                  float* C, const myint* ldc );

int main( int argc, char** argv )
{
    // A is 4x2 embedded in 4x2 array
    // B is 2x3 embedded in 3x3 array
    // C is 4x3 embedded in 5x3 array
    // D = alpha*A*B + beta*C
    myint i, j;
    myint m = 4, n = 3, k = 2, lda = 4, ldb = 3, ldc = 5;
    float alpha = 2, beta = -1;
    float A[ 5*2 ] = { 1, 2, 3, 4,   4, 1, 2, 3 };
    float B[ 3*3 ] = { 1, 3, 0,   2, 1, 0,   3, 2, 0 };
    float C[ 5*3 ] = { 1, 2, 3, 4, 0,   4, 1, 2, 3, 0,   3, 4, 1, 2, 0 };
    float D[ 5*3 ] = { 25, 8, 15, 22, 0,   8, 9, 14, 19, 0,   19, 12, 25, 34 };
    sgemm( "no", "no", &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );
    // check C == D
    for (i = 0; i < ldc*n; ++i) {
        if (C[i] != D[i]) {
            printf( "sgemm failed: C[%d] %.2f != D[%d] %.2f\n",
                    i, C[i], i, D[i] );
            return 1;
        }
    }
    printf( "sgemm ok\n" );
    return 0;
}
