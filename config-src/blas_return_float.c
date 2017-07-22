
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

#define sdot FORTRAN_NAME( sdot, SDOT )

EXTERN_C
float sdot( const myint* n,
            const float* x, const myint* incx,
            const float* y, const myint* incy );

int main( int argc, char** argv )
{
    myint n = 5, ione = 1;
    float x[5] = { 1, 2, 3, 4, 5 };
    float y[5] = { 5, 4, 3, 2, 1 };
    float expect = 35;
    float result = sdot( &n, x, &ione, y, &ione );
    myint okay = (result == expect);
    printf( "sdot result %.2f, expect %.2f, %s\n",
            result, expect, (okay ? "ok" : "failed"));
    return ! okay;
}
